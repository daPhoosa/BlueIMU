#pragma once
// Shared test utilities: double-precision reference math, independent of the library code.

#include "BlueIMU.h"
#include "vendor/doctest.h"
#include <cmath>
#include <new>
#include <initializer_list>
#include <cstring>

constexpr double PI = 3.14159265358979323846;

// Hamilton product a (x) b, computed in double.
inline Quat hamilton(const Quat& a, const Quat& b)
{
   Quat r;
   r.w = (float)((double)a.w*b.w - (double)a.x*b.x - (double)a.y*b.y - (double)a.z*b.z);
   r.x = (float)((double)a.w*b.x + (double)a.x*b.w + (double)a.y*b.z - (double)a.z*b.y);
   r.y = (float)((double)a.w*b.y - (double)a.x*b.z + (double)a.y*b.w + (double)a.z*b.x);
   r.z = (float)((double)a.w*b.z + (double)a.x*b.y - (double)a.y*b.x + (double)a.z*b.w);
   return r;
}

inline Quat quat(double w, double x, double y, double z)
{
   Quat q; q.w = (float)w; q.x = (float)x; q.y = (float)y; q.z = (float)z; return q;
}

// Unit quaternion for a rotation of `angle` rad about unit `axis`.
inline Quat axisAngle(double ax, double ay, double az, double angle)
{
   double n = std::sqrt(ax*ax + ay*ay + az*az), s = std::sin(angle / 2) / n;
   return quat(std::cos(angle / 2), ax*s, ay*s, az*s);
}

// Aerospace ZYX Euler -> Hamilton body-to-world quaternion (double precision).
inline Quat fromEuler(double r, double p, double y)
{
   double cr = std::cos(r/2), sr = std::sin(r/2), cp = std::cos(p/2), sp = std::sin(p/2),
          cy = std::cos(y/2), sy = std::sin(y/2);
   return quat(cr*cp*cy + sr*sp*sy, sr*cp*cy - cr*sp*sy, cr*sp*cy + sr*cp*sy, cr*cp*sy - sr*sp*cy);
}

// Reference rotation of v by unit q (q v q*), via the rotation matrix in double.
inline Vec3 refRotate(const Quat& q, const Vec3& v)
{
   double w=q.w, x=q.x, y=q.y, z=q.z;
   return Vec3((float)((1-2*(y*y+z*z))*v.x + 2*(x*y-w*z)*v.y + 2*(x*z+w*y)*v.z),
               (float)(2*(x*y+w*z)*v.x + (1-2*(x*x+z*z))*v.y + 2*(y*z-w*x)*v.z),
               (float)(2*(x*z-w*y)*v.x + 2*(y*z+w*x)*v.y + (1-2*(x*x+y*y))*v.z));
}

// Angle (rad) between the rotations represented by a and b (sign- and scale-invariant).
inline double rotationDistance(const Quat& a, const Quat& b)
{
   double d = std::fabs((double)a.w*b.w + (double)a.x*b.x + (double)a.y*b.y + (double)a.z*b.z)
            / (std::sqrt((double)a.w*a.w + (double)a.x*a.x + (double)a.y*a.y + (double)a.z*a.z)
             * std::sqrt((double)b.w*b.w + (double)b.x*b.x + (double)b.y*b.y + (double)b.z*b.z));
   return 2.0 * std::acos(d > 1.0 ? 1.0 : d);
}

inline double wrapPi(double a)
{
   while (a >  PI) a -= 2*PI;
   while (a < -PI) a += 2*PI;
   return a;
}

inline double norm(const Quat& q)
{
   return std::sqrt((double)q.w*q.w + (double)q.x*q.x + (double)q.y*q.y + (double)q.z*q.z);
}

// Deterministic pseudo-random doubles (LCG) so failures are reproducible.
struct Rng
{
   uint32_t s;
   explicit Rng(uint32_t seed = 12345u) : s(seed) {}
   double uniform(double lo, double hi) { s = s * 1664525u + 1013904223u; return lo + (hi - lo) * (s >> 8) / 16777216.0; }
   Quat   unitQuat() { return fromEuler(uniform(-PI, PI), uniform(-PI/2, PI/2), uniform(-PI, PI)); }
   Vec3   vec(double m) { return Vec3((float)uniform(-m, m), (float)uniform(-m, m), (float)uniform(-m, m)); }
};

#define CHECK_NEAR(actual, expected, tol)                                                    \
   do {                                                                                      \
      double a_ = (double)(actual), e_ = (double)(expected);                                 \
      INFO(#actual " = " << a_ << ", expected " << e_ << " +/- " << (tol));                  \
      CHECK(std::fabs(a_ - e_) <= (tol));                                                    \
   } while (0)

#define CHECK_VEC_NEAR(a, b, tol)  do { CHECK_NEAR((a).x, (b).x, tol); CHECK_NEAR((a).y, (b).y, tol); CHECK_NEAR((a).z, (b).z, tol); } while (0)

// Quaternions compared component-wise (sign matters).
#define CHECK_QUAT_NEAR(a, b, tol) do { CHECK_NEAR((a).w, (b).w, tol); CHECK_VEC_NEAR(a, b, tol); } while (0)

// ---- simulated time ----
inline void advanceMicros(unsigned long us) { stubMicros() += us; }

// Start both sensor clocks "now" without feeding a sample (so the next input has a valid dt).
inline void primeClocks(BlueIMU& imu)
{
   imu.lastGyroTime  = imu.lastAccelTime  = stubMicros();
   imu.gyroTimeValid = imu.accelTimeValid = true;
}

// Feed one sample dtUs after the previous one.
inline void gyroStep(BlueIMU& imu, unsigned long dtUs, float x, float y, float z)  { advanceMicros(dtUs); imu.inputGyro(x, y, z); }
inline void accelStep(BlueIMU& imu, unsigned long dtUs, float x, float y, float z) { advanceMicros(dtUs); imu.inputAccel(x, y, z); }
inline void gyroStep(BlueIMU& imu, unsigned long dtUs, const Vec3& w)  { gyroStep(imu, dtUs, w.x, w.y, w.z); }
inline void accelStep(BlueIMU& imu, unsigned long dtUs, const Vec3& a) { accelStep(imu, dtUs, a.x, a.y, a.z); }
