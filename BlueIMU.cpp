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
 *	SETTING FUNCTIONS
 */

void BlueIMU::setGyroSampleRate( int rate)
{
   if(rate < 1) rate = 1;
   GyroSamplePeriod  = 1000000UL / rate;  // microseconds per sample
}

void BlueIMU::setMagDeclination( float decDegrees)
{
   float angle = decDegrees * Deg2Rad;
   magneticNorth.x = cos(angle);
   magneticNorth.y = sin(angle);
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

   // Normalize Mag vector
   float magnitude = sqrt(magVectorWorld.x * magVectorWorld.x + magVectorWorld.y * magVectorWorld.y);
   if(magnitude > 0.00001f)
   {
      magVectorWorld = Mul(magVectorWorld, 1.0f / magnitude);

      Vec3 magVectorWorld = Vector(0.0f, 0.0f, atan2(magVectorWorld.y, magVectorWorld.x));

      attitudeEstimate = Mul(Quaternion(magVectorWorld), attitudeEstimate);
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
                               attitudeEstimate.z * attitudeEstimate.z;

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
      velocity = Sum( velocity, Mul(accelVectorWorld, accelDeltaTime) );
      position = Sum( position, Mul(velocity, accelDeltaTime) );
   }

   newAccelData = true;
} 
 
 
void BlueIMU::inputMag(const float& forward, const float& right, const float& down) // Accel input units are in G
{
   // Rotate mag vector from the body frame into the world frame
   Vec3 magVectorWorld = Rotate( attitudeEstimate, Vector(forward, right, down) );

   // Eliminate magnetic inclination
   magVectorWorld.z = 0.0f;

   // Normalize Mag vector
   float magnitudeSquared = magVectorWorld.x * magVectorWorld.x + magVectorWorld.y * magVectorWorld.y;
   if(magnitudeSquared > 0.00001f)
   {
      magVectorWorld = Mul( magVectorWorld, 1.0f / sqrt(magnitudeSquared) );

      // add the mag reading to previous readings (this is effectively a low pass filter)
      // multiple reads are combined into a singe correction
      // Once a drift corection is computed, this will be set to zero
      magAccumulator = Sum( magAccumulator, magVectorWorld );

      newMagData = true;
   }
} 

void velocityCorrection(const float& forward, const float& right, const float& down, const float& weight)
{
   float estimateWeight = 1.0f - weight;
   Vec3 CorrectionVel = Vector(forward, right, down);
   
   velosity = Sum( Mul(velocity, estimateWeight), Mul(CorrectionVel, weight) );
}

void positionCorrection(const float& north,   const float& east,  const float& down, const float& weight)
{
   float estimateWeight = 1.0f - weight;
   Vec3 CorrectionPos = Vector(north, east, down);
   
   velosity = Sum( Mul(position, estimateWeight), Mul(CorrectionPos, weight) );   
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
 
 


 