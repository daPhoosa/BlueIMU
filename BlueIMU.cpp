/*
   BlueIMU  --  Inertial measurement unit library
   Copyright (C) 2016  Phillip J Schmidt

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>
   
 */

#include "BlueIMU.h" 

 
BlueIMU::BlueIMU(int gyroSampleRate)
{
   setGyroSampleRate(gyroSampleRate);
} 


/*
 * SETTING FUNCTIONS
 */

void BlueIMU::setGyroSampleRate( int rate)
{
   if(rate < 1) rate = 1;
   GyroSamplePeriod  = 1000000UL / rate;  // microseconds per sample
}

void BlueIMU::setMagDeclination( float decDegrees)
{
   magDeclination  = decDegrees * Deg2Rad;
   magneticNorth.x = cosf(magDeclination);
   magneticNorth.y = sinf(magDeclination);
}

void BlueIMU::enableVelAndPosEstimation(int accelSampleRate)
{
   if(accelSampleRate < 1) accelSampleRate = 1;
   accelDeltaTime = 1.0f / float(accelSampleRate);
   integrateVelAndPos = true;
}
 
void BlueIMU::setHeadingWithMag(const float& forward,const float& right, const float& down)
{
   // Rotate mag vector from the body frame into the world frame
   Vec3 magVectorWorld = Rotate(attitudeEstimate, Vector(forward, right, down));

   // Eliminate magnetic inclination
   magVectorWorld.z = 0.0f;

   // Only proceed if there is a usable horizontal component
   float magnitudeSquared = magVectorWorld.x * magVectorWorld.x + magVectorWorld.y * magVectorWorld.y;
   if(magnitudeSquared > 1.0e-5f)
   {
      // Current heading of the mag vector in the world frame
      float magHeading = atan2f(magVectorWorld.y, magVectorWorld.x);

      // World-frame yaw rotation that brings the mag vector onto magnetic north (true north + declination)
      Vec3 yawRotation = Vector(0.0f, 0.0f, magDeclination - magHeading);

      // Mul(a, b) is the Hamilton product b*a, so this applies the yaw in the world frame
      attitudeEstimate = Mul(attitudeEstimate, Quaternion(yawRotation));
   }
}
 
 
/*
 *	INPUT FUNCTIONS
 */
 
void BlueIMU::inputGyro(const float& rollRight, const float& pitchUp, const float& yawRight) // Gyro Input units are rad/s
{
	// rotation directions
	// (X) pos roll:  Y --> Z
	// (Y) pos pitch: Z --> X
	// (Z) pos Yaw:   X --> Y
	
	// quaternion integration (rotation composting through multiplication)
	// * The three rotaional rates are converted to a vector structure
	Vec3 gyroVec = Vector(rollRight, pitchUp, yawRight);
	
	if(newDriftCorrection) // check if a drift correction has been generated
   {
		gyroVec = Sum(gyroVec, correctionVectorBody);	// add the correction to the gyro vector to save time doing the correction independently
		newDriftCorrection = false;
	}
	
	// * That vector is combined with the sample time interval to create a quaternion that represents the small attitude change since last update
	rotationDelta = Quaternion(gyroVec, GyroSamplePeriod);

	// * The attitude update quat and the previous quat are combined (through quat multiplication) to create the new attitude estimate
	attitudeEstimate = Mul(rotationDelta, attitudeEstimate);

   normalizeQuat++;
   if(!normalizeQuat)  // re-normalize quaternion magnitude every 256 cycles
   {
      float magnitudeSquared = attitudeEstimate.x * attitudeEstimate.x + 
                               attitudeEstimate.y * attitudeEstimate.y + 
                               attitudeEstimate.z * attitudeEstimate.z +
                               attitudeEstimate.w * attitudeEstimate.w;

      // use fast inverse-square-root since magnitude will always be very close to 1.0
      attitudeEstimate = Mul( attitudeEstimate, (3.0f - magnitudeSquared) * 0.5f ); 
   }
} 
 
 
void BlueIMU::inputAccel(const float& forward, const float& right, const float& down) // Accel input units are in G
{
   // Rotate accel vector from the body frame into the world frame
   Vec3 accelVectorWorld = Rotate( attitudeEstimate, Vector(forward, right, down) );

   // add the accel reading to previous readings (this is effectively a low pass filter)
   // multiple reads are combined into a singe correction
   // Once a correction is computed, this will be set to zero
   accelAccumulator = Sum( accelAccumulator, accelVectorWorld );

   if(integrateVelAndPos)
   {
      // Remove gravity (accelerometer reads VERTICAL at rest) and convert G to m/s^2
      Vec3 linearAccel = Mul( Sum(accelVectorWorld, Mul(VERTICAL, -1.0f)), GRAVITY );

      velocity = Sum( velocity, Mul(linearAccel, accelDeltaTime) );
      position = Sum( position, Mul(velocity, accelDeltaTime) );
   }

   newAccelData = true;
} 
 
 
void BlueIMU::inputMag(const float& forward, const float& right, const float& down) // Mag input units are arbitrary (normalized internally)
{
   // Rotate mag vector from the body frame into the world frame
   Vec3 magVectorWorld = Rotate( attitudeEstimate, Vector(forward, right, down) );

   // Eliminate magnetic inclination
   magVectorWorld.z = 0.0f;

   // Normalize Mag vector
   float magnitudeSquared = magVectorWorld.x * magVectorWorld.x + magVectorWorld.y * magVectorWorld.y;
   if(magnitudeSquared > 1.0e-5f)
   {
      magVectorWorld = Mul( magVectorWorld, 1.0f / sqrtf(magnitudeSquared) );

      // add the mag reading to previous readings (this is effectively a low pass filter)
      // multiple reads are combined into a singe correction
      // Once a drift corection is computed, this will be set to zero
      magAccumulator = Sum( magAccumulator, magVectorWorld );

      newMagData = true;
   }
} 

void BlueIMU::velocityCorrection(const float& north, const float& east, const float& down, const float& weight) // world frame, m/s
{
   float estimateWeight = 1.0f - weight;
   Vec3 CorrectionVel = Vector(north, east, down);
   
   velocity = Sum( Mul(velocity, estimateWeight), Mul(CorrectionVel, weight) );
}

void BlueIMU::positionCorrection(const float& north,   const float& east,  const float& down, const float& weight)
{
   float estimateWeight = 1.0f - weight;
   Vec3 CorrectionPos = Vector(north, east, down);
   
   position = Sum( Mul(position, estimateWeight), Mul(CorrectionPos, weight) );   
}

void BlueIMU::computeDriftCompensation()
{
   Vec3 correctionVectorWorld;

   if(newAccelData)
   {
      // cross product of the accel and the vertical vector creates the pitch/roll error correction vector
      Vec3 pitchRollCorrection = CrossProd(accelAccumulator, VERTICAL);
      correctionVectorWorld.x = pitchRollCorrection.x;
      correctionVectorWorld.y = pitchRollCorrection.y;
      accelAccumulator = ZEROES;  // set accumulator to zero
      newAccelData = false;
   }

   if(newMagData)
   {
      // cross product of the mag and the magnetic north vector creates the yaw error correction vector
      Vec3 yawCorrection = CrossProd(magAccumulator, magneticNorth);
      correctionVectorWorld.z = yawCorrection.z;
      magAccumulator = ZEROES;  // set accumulator to zero
      newMagData = false;
   }

   correctionVectorBody = Rotate(correctionVectorWorld, attitudeEstimate);
   newDriftCorrection = true;
} 


/*
 * OUTPUT FUNCTIONS -- aerospace ZYX Euler angles (NED world, forward/right/down body)
 * roll: right wing down positive, pitch: nose up positive, yaw: clockwise from north positive
 *
 * All angles are taken from rotation-matrix elements written in homogeneous (degree 2) form,
 * e.g. R00 = w^2 + x^2 - y^2 - z^2 rather than 1 - 2(y^2 + z^2). Every atan2 argument pair then
 * scales by the same |q|^2, so the results are exact even if the quaternion has drifted from
 * unit length between renormalizations.
 *
 * Pitch uses atan2(sin, cos) with cos = sqrt(R00^2 + R10^2) instead of asin(), which is
 * ill-conditioned near +/-90 deg (and needs clamping); this form is accurate everywhere.
 */

namespace
{
   // Gimbal lock when cos(pitch) < this (relative to |q|^2): within ~0.057 deg of +/-90 deg.
   // Float rounding in the matrix terms is ~1e-7, so roll/yaw are still accurate (~1e-4 rad)
   // right up to the threshold.
   constexpr float GIMBAL_LOCK_COS = 1.0e-3f;

   inline float fastAtan2(float y, float x)
   {
   #ifdef BLUEIMU_FAST_TRIG
      // Octant reduction + odd minimax polynomial for atan on [0,1]; max error ~1e-5 rad
      float ax = fabsf(x), ay = fabsf(y);
      float mx = (ax > ay) ? ax : ay;
      float mn = (ax > ay) ? ay : ax;
      if(mx == 0.0f) return 0.0f;
      float t = mn / mx;
      float s = t * t;
      float r = t * (0.99997726f + s * (-0.33262347f + s * (0.19354346f + s * (-0.11643287f + s * (0.05265332f + s * -0.01172120f)))));
      if(ay > ax) r = 1.57079637f - r;
      if(x < 0.0f) r = 3.14159274f - r;
      return (y < 0.0f) ? -r : r;
   #else
      return atan2f(y, x);
   #endif
   }

   // Returns true when in gimbal lock. Fills the yaw-row terms (R10, R00) and cos^2(pitch)*|q|^4.
   inline bool yawTerms(const Quat& q, float& r10, float& r00, float& cos2)
   {
      float ww = q.w * q.w, xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
      float n  = ww + xx + yy + zz;
      r10  = 2.0f * (q.x * q.y + q.w * q.z);
      r00  = ww + xx - yy - zz;
      cos2 = r10 * r10 + r00 * r00;
      return cos2 < (GIMBAL_LOCK_COS * GIMBAL_LOCK_COS) * (n * n);
   }

   // Yaw used in gimbal lock (roll forced to 0): atan2(-R01, R11), valid for both +90 and -90 deg
   inline float lockedYaw(const Quat& q)
   {
      return fastAtan2(2.0f * (q.w * q.z - q.x * q.y),
                       q.w * q.w - q.x * q.x + q.y * q.y - q.z * q.z);
   }
}

float BlueIMU::rollRadians() const
{
   const Quat& q = attitudeEstimate;
   float r10, r00, cos2;
   if(yawTerms(q, r10, r00, cos2)) return 0.0f;
   return fastAtan2(2.0f * (q.y * q.z + q.w * q.x),              // R21
                    q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z); // R22
}

float BlueIMU::pitchRadians() const
{
   const Quat& q = attitudeEstimate;
   float r10, r00, cos2;
   yawTerms(q, r10, r00, cos2);
   return fastAtan2(2.0f * (q.w * q.y - q.x * q.z), sqrtf(cos2));  // -R20, cos(pitch)
}

float BlueIMU::yawRadians() const
{
   const Quat& q = attitudeEstimate;
   float r10, r00, cos2;
   if(yawTerms(q, r10, r00, cos2)) return lockedYaw(q);
   return fastAtan2(r10, r00);
}

void BlueIMU::getEulerRadians(float& roll, float& pitch, float& yaw) const
{
   const Quat& q = attitudeEstimate;
   float ww = q.w * q.w, xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
   float wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
   float xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
   float n  = ww + xx + yy + zz;

   float r10  = 2.0f * (xy + wz);
   float r00  = ww + xx - yy - zz;
   float cos2 = r10 * r10 + r00 * r00;

   pitch = fastAtan2(2.0f * (wy - xz), sqrtf(cos2));

   if(cos2 < (GIMBAL_LOCK_COS * GIMBAL_LOCK_COS) * (n * n))
   {
      roll = 0.0f;
      yaw  = fastAtan2(2.0f * (wz - xy), ww - xx + yy - zz);
   }
   else
   {
      roll = fastAtan2(2.0f * (yz + wx), ww - xx - yy + zz);
      yaw  = fastAtan2(r10, r00);
   }
}
 


 