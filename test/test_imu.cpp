// Tests for the BlueIMU class: setup, timing, inputs, corrections and drift compensation

#include "helpers.h"
#include <climits>

TEST_SUITE("setup")
{
   TEST_CASE("constructor initializes state")
   {
      BlueIMU imu;
      CHECK_QUAT_NEAR(imu.attitudeEstimate, quat(1, 0, 0, 0), 0);
      CHECK_VEC_NEAR(imu.velocity, ZEROES, 0);
      CHECK_VEC_NEAR(imu.position, ZEROES, 0);
      CHECK_VEC_NEAR(imu.magneticNorth, Vec3(1, 0, 0), 0);
      CHECK_FALSE(imu.gyroTimeValid);
      CHECK_FALSE(imu.accelTimeValid);
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
      BlueIMU* imu = new (buf) BlueIMU();
      // Inspect raw bytes: reading a garbage bool is UB and may compare as false
      auto rawByte = [](const bool& b) { unsigned char c; std::memcpy(&c, &b, 1); return (int)c; };
      CHECK(rawByte(imu->gyroTimeValid) == 0);
      CHECK(rawByte(imu->accelTimeValid) == 0);
      CHECK(rawByte(imu->newDriftCorrection) == 0);
      CHECK(rawByte(imu->newAccelData) == 0);
      CHECK(rawByte(imu->newMagData) == 0);
      CHECK(rawByte(imu->integrateVelAndPos) == 0);
      CHECK(imu->magDeclination == 0.0f);
      imu->~BlueIMU();
   }

   TEST_CASE("setMagDeclination")
   {
      BlueIMU imu;
      imu.setMagDeclination(90.0f);
      CHECK_NEAR(imu.magDeclination, PI / 2, 1e-6);
      CHECK_VEC_NEAR(imu.magneticNorth, Vec3(0, 1, 0), 1e-6);

      imu.setMagDeclination(-10.0f);
      CHECK_NEAR(imu.magDeclination, -10.0 * PI / 180, 1e-6);
      CHECK_VEC_NEAR(imu.magneticNorth, Vec3((float)std::cos(-10 * PI / 180), (float)std::sin(-10 * PI / 180), 0), 1e-6);
   }

   TEST_CASE("enableVelAndPosEstimation")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      CHECK(imu.integrateVelAndPos);
   }

   TEST_CASE("timing constants")
   {
      CHECK(LARGE_ANGLE_DT_US == 10000UL);
      CHECK(MAX_DT_US == 100000UL);
      CHECK(MAX_DT_US > LARGE_ANGLE_DT_US);
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
         BlueIMU imu;
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
      BlueIMU imu;
      Quat before = fromEuler(0.1, 0.2, 0.3);
      imu.attitudeEstimate = before;
      Vec3 vertical = Rotate(Vec3(0, 0, 1), before);
      imu.setHeadingWithMag(vertical.x, vertical.y, vertical.z);
      CHECK_QUAT_NEAR(imu.attitudeEstimate, before, 0);
   }
}

TEST_SUITE("inputGyro timing")
{
   TEST_CASE("first sample only starts the clock")
   {
      BlueIMU imu;
      imu.correctionVectorBody = Vec3(1, 0, 0);
      imu.newDriftCorrection = true;

      advanceMicros(123456);
      imu.inputGyro(5.0f, 5.0f, 5.0f);
      CHECK(imu.gyroTimeValid);
      CHECK(imu.lastGyroTime == stubMicros());
      CHECK_QUAT_NEAR(imu.attitudeEstimate, quat(1, 0, 0, 0), 0);  // nothing integrated
      CHECK(imu.newDriftCorrection);                                // correction still pending
   }

   TEST_CASE("integrates using the measured interval")
   {
      for (unsigned long dt : {500UL, 2500UL, 10000UL, 20000UL, 80000UL})
      {
         CAPTURE(dt);
         BlueIMU imu;
         primeClocks(imu);
         gyroStep(imu, dt, 0, 0, 2.0f);
         double angle = 2.0 * dt * 1e-6;
         // float rounding + small-angle approximation error (~angle^3/12, only below the threshold)
         CHECK_NEAR(imu.yawRadians(), angle, 1e-5 * angle + angle * angle * angle / 8 + 1e-7);
      }
   }

   TEST_CASE("irregular intervals still integrate the true rotation")
   {
      BlueIMU imu;
      Quat q0 = fromEuler(0.2, 0.5, -1.0);
      imu.attitudeEstimate = q0;
      primeClocks(imu);

      Rng rng(5);
      const double rate = 1.5;   // rad/s about body (1,-2,2)/3
      double total = 0;
      while (total < 2.0)
      {
         // mix of short and long (> LARGE_ANGLE_DT_US) intervals
         unsigned long dt = (rng.uniform(0, 1) < 0.3) ? (unsigned long)rng.uniform(10001, 60000)
                                                     : (unsigned long)rng.uniform(200, 5000);
         gyroStep(imu, dt, (float)(rate / 3), (float)(-2 * rate / 3), (float)(2 * rate / 3));
         total += dt * 1e-6;
      }
      Quat expected = hamilton(q0, axisAngle(1, -2, 2, rate * total));
      CHECK(rotationDistance(imu.attitudeEstimate, expected) < 1e-3);
   }

   TEST_CASE("small-angle quaternion at or below the threshold, large-angle above it")
   {
      Vec3 w(3.0f, -4.0f, 12.0f);   // |w| = 13 rad/s

      BlueIMU imu;
      primeClocks(imu);
      gyroStep(imu, LARGE_ANGLE_DT_US, w);
      CHECK_QUAT_NEAR(imu.rotationDelta, Quaternion(w, LARGE_ANGLE_DT_US), 0);

      BlueIMU imu2;
      primeClocks(imu2);
      gyroStep(imu2, LARGE_ANGLE_DT_US + 1, w);
      CHECK_QUAT_NEAR(imu2.rotationDelta, Quaternion(Mul(w, (LARGE_ANGLE_DT_US + 1) * 1e-6f)), 0);
   }

   TEST_CASE("large intervals are exact (small-angle approximation would not be)")
   {
      // 13 rad/s for 80 ms = 1.04 rad in one step
      Vec3 w(3.0f, -4.0f, 12.0f);
      BlueIMU imu;
      primeClocks(imu);
      gyroStep(imu, 80000UL, w);
      CHECK(rotationDistance(imu.attitudeEstimate, axisAngle(3, -4, 12, 13 * 0.08)) < 1e-5);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.0, 1e-6);

      // Sanity check that this case actually needs the exact form
      CHECK(std::fabs(norm(Quaternion(w, 80000UL)) - 1.0) > 5e-3);
   }

   TEST_CASE("long gaps are clamped to MAX_DT_US")
   {
      BlueIMU imu;
      primeClocks(imu);
      gyroStep(imu, 5000000UL, 0, 0, 1.0f);   // 5 s stall
      CHECK_NEAR(imu.yawRadians(), MAX_DT_US * 1e-6, 1e-5);
   }

   TEST_CASE("handles micros() wrap-around")
   {
      BlueIMU imu;
      stubMicros() = ULONG_MAX - 400;
      imu.inputGyro(0, 0, 0);          // start clock just before the wrap
      gyroStep(imu, 1000, 0, 0, 1.0f); // crosses zero
      CHECK(stubMicros() < 1000);
      CHECK_NEAR(imu.yawRadians(), 0.001, 1e-7);
   }

   TEST_CASE("zero interval changes nothing")
   {
      BlueIMU imu;
      primeClocks(imu);
      imu.inputGyro(5, 5, 5);   // same micros() value
      CHECK_QUAT_NEAR(imu.attitudeEstimate, quat(1, 0, 0, 0), 0);
   }
}

TEST_SUITE("inputGyro")
{
   TEST_CASE("integrates body-frame rates (matches exact rotation)")
   {
      BlueIMU imu;
      Quat q0 = fromEuler(0.2, 0.5, -1.0);
      imu.attitudeEstimate = q0;
      primeClocks(imu);

      const double rate = 1.0;   // rad/s about body (1,-2,2)/3 for 1 s
      for (int i = 0; i < 1000; i++)
         gyroStep(imu, 1000, (float)(rate / 3), (float)(-2 * rate / 3), (float)(2 * rate / 3));

      Quat expected = hamilton(q0, axisAngle(1, -2, 2, rate * 1.0));   // body frame: q0 (x) dq
      CHECK(rotationDistance(imu.attitudeEstimate, expected) < 1e-3);
   }

   TEST_CASE("sign conventions: positive rates give positive Euler angles")
   {
      SUBCASE("roll")  { BlueIMU imu; primeClocks(imu); for (int i = 0; i < 50; i++) gyroStep(imu, 10000, 0.2f, 0, 0); CHECK_NEAR(imu.rollRadians(),  0.1, 1e-4); }
      SUBCASE("pitch") { BlueIMU imu; primeClocks(imu); for (int i = 0; i < 50; i++) gyroStep(imu, 10000, 0, 0.2f, 0); CHECK_NEAR(imu.pitchRadians(), 0.1, 1e-4); }
      SUBCASE("yaw")   { BlueIMU imu; primeClocks(imu); for (int i = 0; i < 50; i++) gyroStep(imu, 10000, 0, 0, 0.2f); CHECK_NEAR(imu.yawRadians(),   0.1, 1e-4); }
   }

   TEST_CASE("zero rate leaves attitude unchanged")
   {
      BlueIMU imu;
      Quat q = fromEuler(0.3, -0.4, 2.0);
      imu.attitudeEstimate = q;
      primeClocks(imu);
      gyroStep(imu, 10000, 0, 0, 0);
      CHECK_QUAT_NEAR(imu.attitudeEstimate, q, 1e-7);
   }

   TEST_CASE("renormalizes the quaternion every 256 samples")
   {
      BlueIMU imu;
      imu.attitudeEstimate = Mul(fromEuler(0.3, 0.2, 0.1), 1.002f);
      primeClocks(imu);
      for (int i = 0; i < 255; i++) gyroStep(imu, 1000, 0, 0, 0);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.002, 1e-5);     // not yet
      gyroStep(imu, 1000, 0, 0, 0);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.0, 1e-5);       // 256th sample
   }

   TEST_CASE("stays normalized over a long run")
   {
      BlueIMU imu;
      primeClocks(imu);
      for (int i = 0; i < 100000; i++) gyroStep(imu, (i % 7) ? 1000 : 15000, 3.0f, -2.0f, 5.0f);
      CHECK_NEAR(norm(imu.attitudeEstimate), 1.0, 1e-4);
   }

   TEST_CASE("applies a pending drift correction exactly once, over the measured interval")
   {
      BlueIMU imu;
      primeClocks(imu);
      imu.correctionVectorBody = Vec3(0.5f, 0, 0);
      imu.newDriftCorrection = true;

      gyroStep(imu, 10000, 0, 0, 0);
      CHECK_FALSE(imu.newDriftCorrection);
      CHECK_NEAR(imu.rollRadians(), 0.005, 1e-6);   // 0.5 rad/s for 10 ms

      gyroStep(imu, 10000, 0, 0, 0);
      CHECK_NEAR(imu.rollRadians(), 0.005, 1e-6);   // not applied again
   }
}

TEST_SUITE("inputAccel")
{
   TEST_CASE("accumulates world-frame readings and flags new data")
   {
      BlueIMU imu;
      Quat q = fromEuler(0.3, -0.2, 1.0);
      imu.attitudeEstimate = q;
      Vec3 a(0.1f, -0.2f, -0.95f);

      accelStep(imu, 1000, a);   // first sample still accumulates for drift correction
      CHECK(imu.newAccelData);
      accelStep(imu, 1000, a);
      CHECK_VEC_NEAR(imu.accelAccumulator, Mul(refRotate(q, a), 2.0f), 1e-5);
   }

   TEST_CASE("does not touch velocity/position when integration is disabled")
   {
      BlueIMU imu;
      for (int i = 0; i < 10; i++) accelStep(imu, 1000, 1, 1, 1);
      CHECK_VEC_NEAR(imu.velocity, ZEROES, 0);
      CHECK_VEC_NEAR(imu.position, ZEROES, 0);
   }

   TEST_CASE("first sample only starts the clock")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      advanceMicros(777777);
      imu.inputAccel(1, 0, -1);
      CHECK(imu.accelTimeValid);
      CHECK_VEC_NEAR(imu.velocity, ZEROES, 0);
      accelStep(imu, 10000, 1, 0, -1);
      CHECK_NEAR(imu.velocity.x, GRAVITY * 0.01, 1e-5);
   }

   TEST_CASE("enabling integration later does not integrate a stale gap")
   {
      BlueIMU imu;
      primeClocks(imu);
      accelStep(imu, 1000, 0, 0, -1);
      advanceMicros(60000000UL);          // a minute passes
      imu.enableVelAndPosEstimation();
      accelStep(imu, 1000, 1, 0, -1);
      // dt is measured from the last sample (clamped), never from the enable call
      CHECK(imu.velocity.x <= GRAVITY * MAX_DT_US * 1e-6 + 1e-4);
   }

   TEST_CASE("removes gravity: a stationary tilted sensor does not move")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      primeClocks(imu);
      Quat q = fromEuler(0.4, 0.3, -2.0);
      imu.attitudeEstimate = q;
      Vec3 g = Rotate(VERTICAL, q);    // what the accelerometer reads at rest
      for (int i = 0; i < 1000; i++) accelStep(imu, 10000, g);
      CHECK(Magnitude(imu.velocity) < 1e-3);
      CHECK(Magnitude(imu.position) < 1e-2);
   }

   TEST_CASE("integrates in m/s and m (discrete semi-implicit Euler)")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      primeClocks(imu);
      const int N = 100; const double dt = 0.01;
      for (int i = 0; i < N; i++) accelStep(imu, 10000, 1, 0, -1);   // level, 1 G forward (north)

      CHECK_NEAR(imu.velocity.x, GRAVITY * N * dt, 1e-3);
      CHECK_NEAR(imu.position.x, GRAVITY * dt * dt * N * (N + 1) / 2.0, 1e-3);
      CHECK_NEAR(imu.velocity.y, 0.0, 1e-6);
      CHECK_NEAR(imu.velocity.z, 0.0, 1e-6);
   }

   TEST_CASE("irregular intervals integrate the measured time")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      primeClocks(imu);
      Rng rng(21);
      double total = 0;
      for (int i = 0; i < 200; i++)
      {
         unsigned long dt = (unsigned long)rng.uniform(500, 30000);
         accelStep(imu, dt, 0.5f, 0, -1);   // 0.5 G north
         total += dt * 1e-6;
      }
      CHECK_NEAR(imu.velocity.x, 0.5 * GRAVITY * total, 1e-3);
   }

   TEST_CASE("long gaps are clamped to MAX_DT_US")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      primeClocks(imu);
      accelStep(imu, 3000000UL, 1, 0, -1);   // 3 s gap
      CHECK_NEAR(imu.velocity.x, GRAVITY * MAX_DT_US * 1e-6, 1e-5);
   }

   TEST_CASE("linear acceleration is resolved into the world frame")
   {
      BlueIMU imu;
      imu.enableVelAndPosEstimation();
      primeClocks(imu);
      Quat q = fromEuler(-0.5, 0.2, 2.5);
      imu.attitudeEstimate = q;
      Vec3 linWorld(1.0f, -2.0f, 0.5f);                                      // m/s^2, NED
      Vec3 specificWorld = Sum(Mul(linWorld, 1.0f / GRAVITY), VERTICAL);    // in G
      Vec3 reading = Rotate(specificWorld, q);
      for (int i = 0; i < 200; i++) accelStep(imu, 5000, reading);   // 1 s
      CHECK_VEC_NEAR(imu.velocity, linWorld, 1e-3);
   }
}

TEST_SUITE("inputMag")
{
   TEST_CASE("accumulates the normalized horizontal component")
   {
      BlueIMU imu;
      imu.inputMag(3, 4, 10);
      CHECK(imu.newMagData);
      CHECK_VEC_NEAR(imu.magAccumulator, Vec3(0.6f, 0.8f, 0), 1e-6);

      imu.inputMag(30, 40, -100);   // magnitude does not matter
      CHECK_VEC_NEAR(imu.magAccumulator, Vec3(1.2f, 1.6f, 0), 1e-6);
   }

   TEST_CASE("is rotated into the world frame")
   {
      BlueIMU imu;
      Quat q = fromEuler(0.2, -0.3, 1.2);
      imu.attitudeEstimate = q;
      Vec3 fieldWorld(0.5f, 0.0f, 0.8f);
      Vec3 m = Rotate(fieldWorld, q);
      imu.inputMag(m.x, m.y, m.z);
      CHECK_VEC_NEAR(imu.magAccumulator, Vec3(1, 0, 0), 1e-5);
   }

   TEST_CASE("ignores a field with no horizontal component")
   {
      BlueIMU imu;
      imu.inputMag(0, 0, 5);
      CHECK_FALSE(imu.newMagData);
      CHECK_VEC_NEAR(imu.magAccumulator, ZEROES, 0);
   }
}

TEST_SUITE("external corrections")
{
   TEST_CASE("velocityCorrection blends in the world frame")
   {
      BlueIMU imu;
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
      BlueIMU imu;
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
      BlueIMU imu;
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
      BlueIMU imu;
      imu.attitudeEstimate = fromEuler(0.3, 0.2, 0.1);
      imu.computeDriftCompensation();
      CHECK_VEC_NEAR(imu.correctionVectorBody, ZEROES, 0);
   }

   TEST_CASE("correct attitude gives a (near) zero correction")
   {
      BlueIMU imu;
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
            BlueIMU imu;
            primeClocks(imu);
            Quat truth = fromEuler(0.1, 0.05, 1.0);
            imu.attitudeEstimate = fromEuler(0.1 + rollErr, 0.05 + pitchErr, 1.0);
            double before = rotationDistance(imu.attitudeEstimate, truth);

            Vec3 g = Rotate(VERTICAL, truth);
            imu.inputAccel(g.x, g.y, g.z);
            imu.computeDriftCompensation();
            gyroStep(imu, 10000, 0, 0, 0);
            CHECK(rotationDistance(imu.attitudeEstimate, truth) < before);
         }
   }

   TEST_CASE("roll and pitch converge to the accelerometer vertical")
   {
      BlueIMU imu;
      primeClocks(imu);
      Quat truth = fromEuler(0.3, -0.4, 0.6);
      imu.attitudeEstimate = fromEuler(-0.2, 0.1, 0.6);   // 0.5 rad roll error, 0.5 rad pitch error
      Vec3 g = Rotate(VERTICAL, truth);

      for (int i = 0; i < 20000; i++)
      {
         gyroStep(imu, 1000, 0, 0, 0);
         accelStep(imu, 0, g);
         if (i % 10 == 9) imu.computeDriftCompensation();
      }
      CHECK_NEAR(imu.rollRadians(),   0.3, 1e-3);
      CHECK_NEAR(imu.pitchRadians(), -0.4, 1e-3);
   }

   TEST_CASE("yaw converges to the magnetometer heading (with declination)")
   {
      BlueIMU imu;
      primeClocks(imu);
      imu.setMagDeclination(8.0f);
      double d = 8.0 * PI / 180;
      Quat truth = fromEuler(0.2, 0.1, 1.0);
      imu.attitudeEstimate = fromEuler(0.2, 0.1, -0.5);   // 1.5 rad yaw error
      Vec3 m = Rotate(Vec3((float)(0.3 * std::cos(d)), (float)(0.3 * std::sin(d)), 0.45f), truth);

      // Correction strength scales with samples accumulated per drift step (no gain term yet),
      // so feed the mag every gyro sample, like the accel test above.
      for (int i = 0; i < 20000; i++)
      {
         gyroStep(imu, 1000, 0, 0, 0);
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
      BlueIMU imu;
      primeClocks(imu);
      Quat truth = fromEuler(0, 0, 0);
      Vec3 g = Rotate(VERTICAL, truth);
      Vec3 m = Rotate(Vec3(0.3f, 0, 0.45f), truth);
      for (int i = 0; i < 20000; i++)
      {
         gyroStep(imu, 1000, 0.01f, -0.01f, 0.01f);   // ~0.57 deg/s bias on every axis
         accelStep(imu, 0, g);
         imu.inputMag(m.x, m.y, m.z);
         if (i % 10 == 9) imu.computeDriftCompensation();
      }
      CHECK(rotationDistance(imu.attitudeEstimate, truth) < 0.05);
   }
}
