// Tests for rollRadians, pitchRadians, yawRadians and getEulerRadians.
// Built and run twice: with libm atan2f and with -DBLUEIMU_FAST_TRIG.

#include "helpers.h"

namespace
{
   const double TOL = 5e-5;   // covers float rounding and the fast-trig approximation (~1e-5)

   void euler(const BlueIMU& imu, double& r, double& p, double& y)
   {
      float rf, pf, yf;
      imu.getEulerRadians(rf, pf, yf);
      r = rf; p = pf; y = yf;
   }
}

TEST_SUITE("Euler outputs")
{
   TEST_CASE("identity attitude is all zeros")
   {
      BlueIMU imu;
      CHECK(imu.rollRadians()  == 0.0f);
      CHECK(imu.pitchRadians() == 0.0f);
      CHECK(imu.yawRadians()   == 0.0f);
   }

   TEST_CASE("single-axis rotations")
   {
      BlueIMU imu;
      for (double a : {-3.0, -1.5, -0.3, 0.3, 1.5, 3.0})
      {
         CAPTURE(a);
         imu.attitudeEstimate = fromEuler(a, 0, 0);
         CHECK_NEAR(imu.rollRadians(), a, TOL);  CHECK_NEAR(imu.pitchRadians(), 0, TOL);  CHECK_NEAR(imu.yawRadians(), 0, TOL);
         imu.attitudeEstimate = fromEuler(0, 0, a);
         CHECK_NEAR(imu.rollRadians(), 0, TOL);  CHECK_NEAR(imu.pitchRadians(), 0, TOL);  CHECK_NEAR(imu.yawRadians(), a, TOL);
         if (std::fabs(a) < PI / 2)
         {
            imu.attitudeEstimate = fromEuler(0, a, 0);
            CHECK_NEAR(imu.rollRadians(), 0, TOL);  CHECK_NEAR(imu.pitchRadians(), a, TOL);  CHECK_NEAR(imu.yawRadians(), 0, TOL);
         }
      }
   }

   TEST_CASE("round-trips random attitudes away from gimbal lock")
   {
      BlueIMU imu;
      Rng rng(7);
      for (int i = 0; i < 20000; i++)
      {
         double r = rng.uniform(-PI, PI), p = rng.uniform(-1.5, 1.5), y = rng.uniform(-PI, PI);
         imu.attitudeEstimate = fromEuler(r, p, y);
         double R, P, Y; euler(imu, R, P, Y);
         INFO("r=" << r << " p=" << p << " y=" << y);
         CHECK_NEAR(wrapPi(R - r), 0, TOL);
         CHECK_NEAR(P - p, 0, TOL);
         CHECK_NEAR(wrapPi(Y - y), 0, TOL);
      }
   }

   TEST_CASE("full quadrant coverage of roll and yaw (atan2 signs)")
   {
      BlueIMU imu;
      for (int k = -11; k <= 12; k++)   // 15 deg steps through (-180, 180]
      {
         double a = k * PI / 12;
         CAPTURE(a);
         imu.attitudeEstimate = fromEuler(a, 0.2, -a);
         CHECK_NEAR(wrapPi(imu.rollRadians() - a), 0, TOL);
         CHECK_NEAR(wrapPi(imu.yawRadians() + a), 0, TOL);
      }
   }

   TEST_CASE("individual functions match getEulerRadians exactly")
   {
      BlueIMU imu;
      Rng rng(11);
      for (int i = 0; i < 5000; i++)
      {
         double p = (i % 2) ? rng.uniform(-PI / 2, PI / 2) : (PI / 2 - std::pow(10.0, -rng.uniform(0, 7))) * (i % 4 ? 1 : -1);
         imu.attitudeEstimate = Mul(fromEuler(rng.uniform(-PI, PI), p, rng.uniform(-PI, PI)), (float)rng.uniform(0.5, 2.0));
         float r, pp, y;
         imu.getEulerRadians(r, pp, y);
         CHECK(imu.rollRadians()  == r);
         CHECK(imu.pitchRadians() == pp);
         CHECK(imu.yawRadians()   == y);
      }
   }

   TEST_CASE("independent of quaternion scale and sign")
   {
      BlueIMU imu;
      Quat q = fromEuler(0.7, -0.6, 2.2);
      for (float s : {0.5f, 0.999f, 1.0f, 1.001f, 3.0f, -1.0f, -2.0f})
      {
         CAPTURE(s);
         imu.attitudeEstimate = Mul(q, s);
         double R, P, Y; euler(imu, R, P, Y);
         CHECK_NEAR(R, 0.7, TOL);
         CHECK_NEAR(P, -0.6, TOL);
         CHECK_NEAR(Y, 2.2, TOL);
      }
   }
}

TEST_SUITE("Euler outputs near gimbal lock")
{
   TEST_CASE("exactly +/-90 deg pitch: roll = 0, combined rotation in yaw")
   {
      BlueIMU imu;
      for (double roll : {-2.0, -0.3, 0.0, 0.3, 1.0})
         for (double yaw : {-2.5, 0.0, 1.0})
         {
            CAPTURE(roll); CAPTURE(yaw);
            imu.attitudeEstimate = fromEuler(roll, PI / 2, yaw);
            double R, P, Y; euler(imu, R, P, Y);
            CHECK(R == 0.0);
            CHECK_NEAR(P, PI / 2, 1e-6);
            CHECK_NEAR(wrapPi(Y - (yaw - roll)), 0, TOL);   // nose up: yaw - roll

            imu.attitudeEstimate = fromEuler(roll, -PI / 2, yaw);
            euler(imu, R, P, Y);
            CHECK(R == 0.0);
            CHECK_NEAR(P, -PI / 2, 1e-6);
            CHECK_NEAR(wrapPi(Y - (yaw + roll)), 0, TOL);   // nose down: yaw + roll
         }
   }

   TEST_CASE("reported angles always describe the true rotation")
   {
      BlueIMU imu;
      Rng rng(3);
      for (int i = 0; i < 20000; i++)
      {
         double sign = (i % 2) ? 1.0 : -1.0;
         double p = sign * (PI / 2 - std::pow(10.0, -rng.uniform(0, 8)));
         Quat q = fromEuler(rng.uniform(-PI, PI), p, rng.uniform(-PI, PI));
         imu.attitudeEstimate = q;
         double R, P, Y; euler(imu, R, P, Y);
         REQUIRE(std::isfinite(R)); REQUIRE(std::isfinite(P)); REQUIRE(std::isfinite(Y));
         CHECK(std::fabs(P) <= PI / 2 + 1e-6);
         // inside the lock band roll is forced to 0, costing at most ~2x the band width (~2e-3 rad)
         CHECK(rotationDistance(fromEuler(R, P, Y), q) < 3e-3);
      }
   }

   TEST_CASE("pitch is monotonic and accurate approaching 90 deg")
   {
      BlueIMU imu;
      double last = -1;
      for (int k = 0; k <= 2000; k++)
      {
         double p = 1.4 + (PI / 2 - 1.4) * k / 2000.0;
         imu.attitudeEstimate = fromEuler(0.5, p, -1.0);
         double P = imu.pitchRadians();
         CHECK(P >= last);
         CHECK_NEAR(P, p, 1e-5);
         last = P;
      }
   }

   TEST_CASE("degenerate zero quaternion does not produce NaN")
   {
      BlueIMU imu;
      imu.attitudeEstimate = quat(0, 0, 0, 0);
      double R, P, Y; euler(imu, R, P, Y);
      CHECK(std::isfinite(R)); CHECK(std::isfinite(P)); CHECK(std::isfinite(Y));
      CHECK(std::isfinite(imu.rollRadians()));
      CHECK(std::isfinite(imu.pitchRadians()));
      CHECK(std::isfinite(imu.yawRadians()));
   }
}
