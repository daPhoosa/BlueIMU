// Tests for the inline math helpers and constants in BlueIMU.h

#include "helpers.h"

TEST_SUITE("constants")
{
   TEST_CASE("angle conversions")
   {
      CHECK_NEAR(Deg2Rad * 180.0f, PI, 1e-6);
      CHECK_NEAR(Rad2Deg * PI, 180.0, 1e-4);
      CHECK_NEAR(Deg2Rad * Rad2Deg, 1.0, 1e-6);
   }

   TEST_CASE("GRAVITY, VERTICAL, ZEROES")
   {
      CHECK(GRAVITY == doctest::Approx(9.80665));
      CHECK_VEC_NEAR(VERTICAL, Vec3(0, 0, -1), 0);
      CHECK_VEC_NEAR(ZEROES, Vec3(0, 0, 0), 0);
   }
}

TEST_SUITE("Vec3 / Quat construction")
{
   TEST_CASE("default values")
   {
      Vec3 v;
      CHECK_VEC_NEAR(v, Vec3(0, 0, 0), 0);
      Quat q;
      CHECK_QUAT_NEAR(q, quat(1, 0, 0, 0), 0);
   }

   TEST_CASE("Vector() and Vec3 constructor")
   {
      Vec3 v = Vector(1.5f, -2.0f, 3.25f);
      CHECK(v.x == 1.5f);
      CHECK(v.y == -2.0f);
      CHECK(v.z == 3.25f);
      CHECK_VEC_NEAR(Vec3(1.5f, -2.0f, 3.25f), v, 0);
   }

   TEST_CASE("Vector() is usable in constant expressions")
   {
      constexpr Vec3 v = Vector(1.0f, 2.0f, 3.0f);
      static_assert(v.y == 2.0f, "constexpr Vector");
      CHECK(v.z == 3.0f);
   }
}

TEST_SUITE("vector arithmetic")
{
   TEST_CASE("Sum")
   {
      CHECK_VEC_NEAR(Sum(Vec3(1, 2, 3), Vec3(-4, 0.5f, 10)), Vec3(-3, 2.5f, 13), 0);
   }

   TEST_CASE("Mul(Vec3, float) and Mul(float, Vec3)")
   {
      CHECK_VEC_NEAR(Mul(Vec3(1, -2, 3), 2.0f), Vec3(2, -4, 6), 0);
      CHECK_VEC_NEAR(Mul(-0.5f, Vec3(1, -2, 3)), Vec3(-0.5f, 1, -1.5f), 0);
   }

   TEST_CASE("Magnitude(Vec3)")
   {
      CHECK_NEAR(Magnitude(Vec3(3, 4, 12)), 13.0, 1e-6);
      CHECK(Magnitude(Vec3(0, 0, 0)) == 0.0f);
   }

   TEST_CASE("CrossProd")
   {
      CHECK_VEC_NEAR(CrossProd(Vec3(1, 0, 0), Vec3(0, 1, 0)), Vec3(0, 0, 1), 0);
      CHECK_VEC_NEAR(CrossProd(Vec3(0, 1, 0), Vec3(0, 0, 1)), Vec3(1, 0, 0), 0);
      CHECK_VEC_NEAR(CrossProd(Vec3(0, 0, 1), Vec3(1, 0, 0)), Vec3(0, 1, 0), 0);

      Rng rng;
      for (int i = 0; i < 100; i++)
      {
         Vec3 a = rng.vec(5), b = rng.vec(5);
         Vec3 c = CrossProd(a, b);
         // anti-commutative and orthogonal to both inputs
         CHECK_VEC_NEAR(CrossProd(b, a), Mul(c, -1.0f), 1e-5);
         CHECK_NEAR(c.x*a.x + c.y*a.y + c.z*a.z, 0.0, 1e-4);
         CHECK_NEAR(c.x*b.x + c.y*b.y + c.z*b.z, 0.0, 1e-4);
      }
   }
}

TEST_SUITE("quaternion arithmetic")
{
   TEST_CASE("Mul(Quat, Quat) is the Hamilton product with operands swapped")
   {
      // Library convention: Mul(lhs, rhs) == rhs (x) lhs
      Rng rng;
      for (int i = 0; i < 100; i++)
      {
         Quat a = rng.unitQuat(), b = rng.unitQuat();
         CHECK_QUAT_NEAR(Mul(a, b), hamilton(b, a), 1e-6);
      }
   }

   TEST_CASE("Mul(Quat, Quat) identity and composition order")
   {
      Quat id;
      Quat a = axisAngle(1, 2, 3, 0.7);
      CHECK_QUAT_NEAR(Mul(a, id), a, 0);
      CHECK_QUAT_NEAR(Mul(id, a), a, 0);

      // Mul(a, b) applies a first, then b
      Quat b = axisAngle(-1, 0, 2, 1.3);
      Vec3 v(0.3f, -1.2f, 2.0f);
      CHECK_VEC_NEAR(Rotate(Mul(a, b), v), Rotate(b, Rotate(a, v)), 1e-5);
   }

   TEST_CASE("Mul(Quat, Vec3) treats the vector as a pure quaternion")
   {
      Rng rng;
      for (int i = 0; i < 50; i++)
      {
         Quat q = rng.unitQuat();
         Vec3 v = rng.vec(3);
         CHECK_QUAT_NEAR(Mul(q, v), Mul(q, quat(0, v.x, v.y, v.z)), 1e-6);
      }
   }

   TEST_CASE("Mul(Quat, float)")
   {
      CHECK_QUAT_NEAR(Mul(quat(1, -2, 3, -4), 0.5f), quat(0.5, -1, 1.5, -2), 0);
   }

   TEST_CASE("Magnitude(Quat)")
   {
      CHECK_NEAR(Magnitude(quat(1, 2, 2, 4)), 5.0, 1e-6);
      CHECK_NEAR(Magnitude(axisAngle(1, 1, 1, 2.0)), 1.0, 1e-6);
   }
}

TEST_SUITE("rotation quaternions")
{
   TEST_CASE("Quaternion(rotationVector): exact axis-angle")
   {
      Quat q = Quaternion(Vec3(0, 0, (float)(PI / 2)));
      CHECK_QUAT_NEAR(q, axisAngle(0, 0, 1, PI / 2), 1e-6);
      CHECK_VEC_NEAR(Rotate(q, Vec3(1, 0, 0)), Vec3(0, 1, 0), 1e-6);  // +yaw: X -> Y

      // large angle about an arbitrary axis
      Quat big = Quaternion(Mul(Vec3(0.6f, 0.0f, 0.8f), (float)(0.9 * PI)));
      CHECK_QUAT_NEAR(big, axisAngle(0.6, 0, 0.8, 0.9 * PI), 1e-6);
      CHECK_NEAR(Magnitude(big), 1.0, 1e-6);
   }

   TEST_CASE("Quaternion(rotationVector): negligible rotation returns identity")
   {
      CHECK_QUAT_NEAR(Quaternion(Vec3(0, 0, 0)), quat(1, 0, 0, 0), 0);
      CHECK_QUAT_NEAR(Quaternion(Vec3(1e-8f, 0, 0)), quat(1, 0, 0, 0), 0);
   }

   TEST_CASE("Quaternion(rate, microseconds): small-angle approximation")
   {
      // 1 rad/s about (1,2,2)/3 for 1 ms -> 0.001 rad rotation
      Vec3 rate(1.0f / 3, 2.0f / 3, 2.0f / 3);
      Quat q = Quaternion(rate, 1000UL);
      CHECK_QUAT_NEAR(q, axisAngle(1, 2, 2, 0.001), 1e-7);
      CHECK_NEAR(Magnitude(q), 1.0, 1e-6);

      // zero time or zero rate is the identity
      CHECK_QUAT_NEAR(Quaternion(rate, 0UL), quat(1, 0, 0, 0), 0);
      CHECK_QUAT_NEAR(Quaternion(Vec3(0, 0, 0), 5000UL), quat(1, 0, 0, 0), 0);
   }
}

TEST_SUITE("vector rotation")
{
   TEST_CASE("Rotate(q, v) matches the rotation matrix (body -> world)")
   {
      Rng rng;
      for (int i = 0; i < 200; i++)
      {
         Quat q = rng.unitQuat();
         Vec3 v = rng.vec(10);
         CHECK_VEC_NEAR(Rotate(q, v), refRotate(q, v), 1e-5);
      }
   }

   TEST_CASE("Rotate(v, q) is the inverse rotation (world -> body)")
   {
      Rng rng(99);
      for (int i = 0; i < 200; i++)
      {
         Quat q = rng.unitQuat();
         Vec3 v = rng.vec(10);
         CHECK_VEC_NEAR(Rotate(Rotate(q, v), q), v, 2e-5);
         CHECK_VEC_NEAR(Rotate(q, Rotate(v, q)), v, 2e-5);
         CHECK_NEAR(Magnitude(Rotate(q, v)), Magnitude(v), 1e-5);
      }
   }

   TEST_CASE("Rotate: sign conventions for roll, pitch, yaw")
   {
      // +roll (right wing down): body right axis points down in the world
      CHECK_VEC_NEAR(Rotate(fromEuler(PI / 2, 0, 0), Vec3(0, 1, 0)), Vec3(0, 0, 1), 1e-6);
      // +pitch (nose up): body forward axis points up (-Z in NED)
      CHECK_VEC_NEAR(Rotate(fromEuler(0, PI / 2, 0), Vec3(1, 0, 0)), Vec3(0, 0, -1), 1e-6);
      // +yaw (nose right): body forward axis points east
      CHECK_VEC_NEAR(Rotate(fromEuler(0, 0, PI / 2), Vec3(1, 0, 0)), Vec3(0, 1, 0), 1e-6);
   }
}
