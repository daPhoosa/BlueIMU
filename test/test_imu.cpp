// Tests for the BlueIMU class: setup, inputs, corrections and drift compensation

#include "helpers.h"

TEST_SUITE("setup")
{
   TEST_CASE("constructor initializes state")
   {
      BlueIMU imu(500);
      CHECK(imu.GyroSamplePeriod == 2000UL);
      CHECK_QUAT_NEAR(imu.attitudeEstimate, quat(1, 0, 0, 0), 0);
      CHECK_VEC_NEAR(imu.velocity, ZEROES, 0);
      CHECK_VEC_NEAR(imu.position, ZEROES, 0);
      CHECK_VEC_NEAR(imu.magneticNorth, Vec3(1, 0, 0), 0);
      CHECK_FALSE(imu.newDriftCorrection);
      CHECK_FALSE(imu.newAccelData);
      CHECK_FALSE(imu.newMagData);
      CHECK_FALSE(imu.integrateVelAndPos);
      CHECK(imu.magDeclination == 0.0f);
   }

   TEST_CASE("constructor initializes flags even in dirty memory")
   {
      // Placement-new over garbage: members must not inherit whatever was there before
      alignas(BlueIMU) unsigned char buf[sizeof(BlueIMU)];
      for (unsigned i = 0; i < sizeof(buf); i++) buf[i] = 0xA5;
      BlueIMU* imu = new (buf) BlueIMU(100);
      // Inspect raw bytes: reading a garbage bool is UB and may compare as false
      auto rawByte = [](const bool& b) { unsigned char c; std::memcpy(&c, &b, 1); return (int)c; };
      CHECK(rawByte(imu->newDriftCorrection) == 0);
      CHECK(rawByte(imu->newAccelData) == 0);
      CHECK(rawByte(imu->newMagData) == 0);
      CHECK(rawByte(imu->integrateVelAndPos) == 0);
      CHECK(imu->magDeclination == 0.0f);
      CHECK(imu->accelDeltaTime == 0.0f);
      imu->~BlueIMU();
   }

   TEST_CASE("setGyroSampleRate")
   {
      BlueIMU imu(100);
      CHECK(imu.GyroSamplePeriod == 10000UL);
      imu.setGyroSampleRate(1000);
      CHECK(imu.GyroSamplePeriod == 1000UL);
      imu.setGyroSampleRate(0);    // clamped to 1 Hz
      CHECK(imu.GyroSamplePeriod == 1000000UL);
      imu.setGyroSampleRate(-5);
      CHECK(imu.GyroSamplePeriod == 1000000UL);
   }

   TEST_CASE("setMagDeclination")
   {
      BlueIMU imu(100);
      imu.setMagDeclination(90.0f);
      CHECK_NEAR(imu.magDeclination, PI / 2, 1e-6);
      CHECK_VEC_NEAR(imu.magneticNorth, Vec3(0, 1, 0), 1e-6);

      imu.setMagDeclination(-10.0f);
      CHECK_NEAR(imu.magDeclination, -10.0 * PI / 180, 1e-6);
      CHECK_VEC_NEAR(imu.magneticNorth, Vec3((float)std::cos(-10 * PI / 180), (float)std::sin(-10 * PI / 180), 0), 1e-6);
   }

   TEST_CASE("enableVelAndPosEstimation")
   {
      BlueIMU imu(100);
      imu.enableVelAndPosEstimation(200);
      CHECK(imu.integrateVelAndPos);
      CHECK_NEAR(imu.accelDeltaTime, 0.005, 1e-9);
      imu.enableVelAndPosEstimation(0);   // clamped to 1 Hz
      CHECK_NEAR(imu.accelDeltaTime, 1.0, 0);
   }
}

TEST_SUITE("setHeadingWithMag")
{
   // Magnetometer reading in the body frame for a given true attitude and declination.
   static Vec3 magReading(const Quat& truth, double declinationDeg)
   {
      double d = declinationDeg * PI / 180;
      Vec3 fieldWorld((float)(0.3 * std::cos(d)), (float)(0.3 * std::sin(d)), 0.45f);  // inclined field
      return Rotate(fieldWorld, truth);
   }

   TEST_CASE("aligns heading to true north, keeps tilt, applies declination")
   {
      for (double decl : {0.0, 12.5, -20.0})
      {
         CAPTURE(decl);
         BlueIMU imu(100);
         imu.setMagDeclination((float)decl);
         Quat truth = fromEuler(0.4, -0.25, 0.7);
         imu.attitudeEstimate = fromEuler(0.4, -0.25, -2.0);   // right tilt, wrong heading

         Vec3 m = magReading(truth, decl);
         imu.setHeadingWithMag(m.x, m.y, m.z);

         CHECK(rotationDistance(imu.attitudeEstimate, truth) < 1e-4);
         CHECK_NEAR(imu.yawRadians(),   0.7,  1e-4);
         CHECK_NEAR(imu.rollRadians(),  0.4,  1e-4);
         CHECK_NEAR(imu.pitchRadians(), -0.25, 1e-4);
      }
   }

   TEST_CASE("field with no horizontal component leaves attitude unchanged")
   {
      BlueIMU imu(100);
      Quat before = fromEuler(0.1, 0.2, 0.3);
      imu.attitudeEstimate = before;
      Vec3 vertical = Rotate(Vec3(0, 0, 1), before);
      imu.setHeadingWithMag(vertical.x, vertical.y, vertical.z);
      CHECK_QUAT_NEAR(imu.attitudeEstimate, before, 0);
   }
}

TEST_SUITE("inputGyro")
{
   TEST_CASE("integrates body-frame rates (matches exact rotation)")
   {
      BlueIMU imu(1000);
      Quat q0 = fromEuler(0.2, 0.5, -1.0);
      imu.attitudeEstimate = q0;

      const double rate = 1.0;   // rad/s about body (1,-2,2)/3 for 1 s
      for (int i = 0; i < 1000; i++)
         imu.inputGyro((float)(rate / 3), (float)(-2 * rate / 3), (float)(2 * rate / 3));

      Quat expected = hamilton(q0, axisAngle(1, -2, 2, rate * 1.0));   // body frame: q0 (x) dq
      CHECK(rotationDistance(imu.attitudeEstimate, expected) < 1e-3);
   }

   TEST_CASE("sign conventions: positive rates give positive Euler angles")
   {
      SUBCASE("roll")  { BlueIMU imu(100); for (int i = 0; i < 50; i++) imu.inputGyro(0.2f, 0, 0); CHECK_NEAR(imu.rollRadians(),  0.1, 1e-4); }
      SUBCASE("pitch") { BlueIMU imu(100); for (int i = 0; i < 50; i++) imu.inputGyro(0, 0.2f, 0); CHECK_NEAR(imu.pitchRadians(), 0.1, 1e-4); }
      SUBCASE("yaw")   { BlueIMU imu(100); for (int i = 0; i < 50; i++) imu.inputGyro(0, 0, 0.2f); CHECK_NEAR(imu.yawRadians(),   0.1, 1e-4); }
   }

   TEST_CASE("zero rate leaves attitude unchanged")
   {
      BlueIMU imu(100);
      Quat q = fromEuler(0.3, -0.4, 2.0);
      imu.attitudeEstimate = q;
      imu.inputGyro(0, 0, 0);
      CHECK_QUAT_NEAR(imu.attitudeEstimate, q, 1e-7);
   }

   TEST_CASE("renormalizes the quaternion every 256 samples")
   {
      BlueIMU imu(100);
      imu.attitudeEstimate = Mul(fromEuler(0.3, 0.2, 0.1), 1.002f);
      for (int i = 0; i < 255; i++) imu.inputGyro(0, 0, 0);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.002, 1e-5);     // not yet
      imu.inputGyro(0, 0, 0);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.0, 1e-5);       // 256th sample
   }

   TEST_CASE("stays normalized over a long run")
   {
      BlueIMU imu(1000);
      for (int i = 0; i < 100000; i++) imu.inputGyro(3.0f, -2.0f, 5.0f);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.0, 1e-4);
   }

   TEST_CASE("applies a pending drift correction exactly once")
   {
      BlueIMU imu(100);
      imu.correctionVectorBody = Vec3(0.5f, 0, 0);
      imu.newDriftCorrection = true;

      imu.inputGyro(0, 0, 0);
      CHECK_FALSE(imu.newDriftCorrection);
      CHECK_NEAR(imu.rollRadians(), 0.005, 1e-6);   // 0.5 rad/s for 10 ms

      imu.inputGyro(0, 0, 0);
      CHECK_NEAR(imu.rollRadians(), 0.005, 1e-6);   // not applied again
   }
}

TEST_SUITE("inputAccel")
{
   TEST_CASE("accumulates world-frame readings and flags new data")
   {
      BlueIMU imu(100);
      Quat q = fromEuler(0.3, -0.2, 1.0);
      imu.attitudeEstimate = q;
      Vec3 a(0.1f, -0.2f, -0.95f);

      imu.inputAccel(a.x, a.y, a.z);
      CHECK(imu.newAccelData);
      imu.inputAccel(a.x, a.y, a.z);
      CHECK_VEC_NEAR(imu.accelAccumulator, Mul(refRotate(q, a), 2.0f), 1e-5);
   }

   TEST_CASE("does not touch velocity/position when integration is disabled")
   {
      BlueIMU imu(100);
      for (int i = 0; i < 10; i++) imu.inputAccel(1, 1, 1);
      CHECK_VEC_NEAR(imu.velocity, ZEROES, 0);
      CHECK_VEC_NEAR(imu.position, ZEROES, 0);
   }

   TEST_CASE("removes gravity: a stationary tilted sensor does not move")
   {
      BlueIMU imu(100);
      imu.enableVelAndPosEstimation(100);
      Quat q = fromEuler(0.4, 0.3, -2.0);
      imu.attitudeEstimate = q;
      Vec3 g = Rotate(VERTICAL, q);    // what the accelerometer reads at rest
      for (int i = 0; i < 1000; i++) imu.inputAccel(g.x, g.y, g.z);
      CHECK(Magnitude(imu.velocity) < 1e-3);
      CHECK(Magnitude(imu.position) < 1e-2);
   }

   TEST_CASE("integrates in m/s and m (discrete semi-implicit Euler)")
   {
      BlueIMU imu(100);
      imu.enableVelAndPosEstimation(100);
      const int N = 100; const double dt = 0.01;
      for (int i = 0; i < N; i++) imu.inputAccel(1, 0, -1);   // level, 1 G forward (north)

      CHECK_NEAR(imu.velocity.x, GRAVITY * N * dt, 1e-3);
      CHECK_NEAR(imu.position.x, GRAVITY * dt * dt * N * (N + 1) / 2.0, 1e-3);
      CHECK_NEAR(imu.velocity.y, 0.0, 1e-6);
      CHECK_NEAR(imu.velocity.z, 0.0, 1e-6);
   }

   TEST_CASE("linear acceleration is resolved into the world frame")
   {
      BlueIMU imu(200);
      imu.enableVelAndPosEstimation(200);
      Quat q = fromEuler(-0.5, 0.2, 2.5);
      imu.attitudeEstimate = q;
      Vec3 linWorld(1.0f, -2.0f, 0.5f);                                      // m/s^2, NED
      Vec3 specificWorld = Sum(Mul(linWorld, 1.0f / GRAVITY), VERTICAL);    // in G
      Vec3 reading = Rotate(specificWorld, q);
      for (int i = 0; i < 200; i++) imu.inputAccel(reading.x, reading.y, reading.z);   // 1 s
      CHECK_VEC_NEAR(imu.velocity, linWorld, 1e-3);
   }
}

TEST_SUITE("inputMag")
{
   TEST_CASE("accumulates the normalized horizontal component")
   {
      BlueIMU imu(100);
      imu.inputMag(3, 4, 10);
      CHECK(imu.newMagData);
      CHECK_VEC_NEAR(imu.magAccumulator, Vec3(0.6f, 0.8f, 0), 1e-6);

      imu.inputMag(30, 40, -100);   // magnitude does not matter
      CHECK_VEC_NEAR(imu.magAccumulator, Vec3(1.2f, 1.6f, 0), 1e-6);
   }

   TEST_CASE("is rotated into the world frame")
   {
      BlueIMU imu(100);
      Quat q = fromEuler(0.2, -0.3, 1.2);
      imu.attitudeEstimate = q;
      Vec3 fieldWorld(0.5f, 0.0f, 0.8f);
      Vec3 m = Rotate(fieldWorld, q);
      imu.inputMag(m.x, m.y, m.z);
      CHECK_VEC_NEAR(imu.magAccumulator, Vec3(1, 0, 0), 1e-5);
   }

   TEST_CASE("ignores a field with no horizontal component")
   {
      BlueIMU imu(100);
      imu.inputMag(0, 0, 5);
      CHECK_FALSE(imu.newMagData);
      CHECK_VEC_NEAR(imu.magAccumulator, ZEROES, 0);
   }
}

TEST_SUITE("external corrections")
{
   TEST_CASE("velocityCorrection blends in the world frame")
   {
      BlueIMU imu(100);
      imu.velocity = Vec3(10, -4, 2);
      imu.attitudeEstimate = fromEuler(0, 0, PI / 2);   // attitude must not affect the blend

      imu.velocityCorrection(2, 0, 6, 0.0f);
      CHECK_VEC_NEAR(imu.velocity, Vec3(10, -4, 2), 0);
      imu.velocityCorrection(2, 0, 6, 0.25f);
      CHECK_VEC_NEAR(imu.velocity, Vec3(8, -3, 3), 1e-5);
      imu.velocityCorrection(2, 0, 6, 1.0f);
      CHECK_VEC_NEAR(imu.velocity, Vec3(2, 0, 6), 1e-6);
   }

   TEST_CASE("positionCorrection blends in the world frame")
   {
      BlueIMU imu(100);
      imu.position = Vec3(100, 50, -10);

      imu.positionCorrection(0, 0, 0, 0.0f);
      CHECK_VEC_NEAR(imu.position, Vec3(100, 50, -10), 0);
      imu.positionCorrection(0, 0, 0, 0.5f);
      CHECK_VEC_NEAR(imu.position, Vec3(50, 25, -5), 1e-5);
      imu.positionCorrection(-1, 2, 3, 1.0f);
      CHECK_VEC_NEAR(imu.position, Vec3(-1, 2, 3), 1e-6);
   }
}

TEST_SUITE("computeDriftCompensation")
{
   TEST_CASE("clears accumulators and flags, sets a pending correction")
   {
      BlueIMU imu(100);
      imu.inputAccel(0.1f, 0, -1);
      imu.inputMag(1, 0.2f, 0.5f);
      imu.computeDriftCompensation();
      CHECK(imu.newDriftCorrection);
      CHECK_FALSE(imu.newAccelData);
      CHECK_FALSE(imu.newMagData);
      CHECK_VEC_NEAR(imu.accelAccumulator, ZEROES, 0);
      CHECK_VEC_NEAR(imu.magAccumulator, ZEROES, 0);
   }

   TEST_CASE("no new data gives a zero correction")
   {
      BlueIMU imu(100);
      imu.attitudeEstimate = fromEuler(0.3, 0.2, 0.1);
      imu.computeDriftCompensation();
      CHECK_VEC_NEAR(imu.correctionVectorBody, ZEROES, 0);
   }

   TEST_CASE("correct attitude gives a (near) zero correction")
   {
      BlueIMU imu(100);
      Quat q = fromEuler(0.3, -0.2, 0.8);
      imu.attitudeEstimate = q;
      Vec3 g = Rotate(VERTICAL, q);
      Vec3 m = Rotate(Vec3(0.4f, 0, 0.6f), q);
      imu.inputAccel(g.x, g.y, g.z);
      imu.inputMag(m.x, m.y, m.z);
      imu.computeDriftCompensation();
      CHECK(Magnitude(imu.correctionVectorBody) < 1e-5);
   }

   TEST_CASE("a single correction step reduces roll/pitch error")
   {
      for (double rollErr : {0.2, -0.2})
         for (double pitchErr : {0.15, -0.15})
         {
            CAPTURE(rollErr); CAPTURE(pitchErr);
            BlueIMU imu(100);
            Quat truth = fromEuler(0.1, 0.05, 1.0);
            imu.attitudeEstimate = fromEuler(0.1 + rollErr, 0.05 + pitchErr, 1.0);
            double before = rotationDistance(imu.attitudeEstimate, truth);

            Vec3 g = Rotate(VERTICAL, truth);
            imu.inputAccel(g.x, g.y, g.z);
            imu.computeDriftCompensation();
            imu.inputGyro(0, 0, 0);
            CHECK(rotationDistance(imu.attitudeEstimate, truth) < before);
         }
   }

   TEST_CASE("roll and pitch converge to the accelerometer vertical")
   {
      BlueIMU imu(1000);
      Quat truth = fromEuler(0.3, -0.4, 0.6);
      imu.attitudeEstimate = fromEuler(-0.2, 0.1, 0.6);   // 0.5 rad roll error, 0.5 rad pitch error
      Vec3 g = Rotate(VERTICAL, truth);

      for (int i = 0; i < 20000; i++)
      {
         imu.inputGyro(0, 0, 0);
         imu.inputAccel(g.x, g.y, g.z);
         if (i % 10 == 9) imu.computeDriftCompensation();
      }
      CHECK_NEAR(imu.rollRadians(),   0.3, 1e-3);
      CHECK_NEAR(imu.pitchRadians(), -0.4, 1e-3);
   }

   TEST_CASE("yaw converges to the magnetometer heading (with declination)")
   {
      BlueIMU imu(1000);
      imu.setMagDeclination(8.0f);
      double d = 8.0 * PI / 180;
      Quat truth = fromEuler(0.2, 0.1, 1.0);
      imu.attitudeEstimate = fromEuler(0.2, 0.1, -0.5);   // 1.5 rad yaw error
      Vec3 m = Rotate(Vec3((float)(0.3 * std::cos(d)), (float)(0.3 * std::sin(d)), 0.45f), truth);

      // Correction strength scales with samples accumulated per drift step (no gain term yet),
      // so feed the mag every gyro sample, like the accel test above.
      for (int i = 0; i < 20000; i++)
      {
         imu.inputGyro(0, 0, 0);
         imu.inputMag(m.x, m.y, m.z);
         if (i % 10 == 9) imu.computeDriftCompensation();
      }
      CHECK_NEAR(wrapPi(imu.yawRadians() - 1.0), 0.0, 1e-3);
      CHECK_NEAR(imu.rollRadians(),  0.2, 1e-3);
      CHECK_NEAR(imu.pitchRadians(), 0.1, 1e-3);
   }

   TEST_CASE("cancels a constant gyro bias to a bounded error")
   {
      // No integral term, so a small steady-state offset remains; it must stay bounded.
      BlueIMU imu(1000);
      Quat truth = fromEuler(0, 0, 0);
      Vec3 g = Rotate(VERTICAL, truth);
      Vec3 m = Rotate(Vec3(0.3f, 0, 0.45f), truth);
      for (int i = 0; i < 20000; i++)
      {
         imu.inputGyro(0.01f, -0.01f, 0.01f);   // ~0.57 deg/s bias on every axis
         imu.inputAccel(g.x, g.y, g.z);
         imu.inputMag(m.x, m.y, m.z);
         if (i % 10 == 9) imu.computeDriftCompensation();
      }
      CHECK(rotationDistance(imu.attitudeEstimate, truth) < 0.05);
   }
}
