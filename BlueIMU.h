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



#ifndef BlueIMU_h
	#define BlueIMU_h
	
	#include <arduino.h>
	
   struct Quat
   {
      float w = 1.0f;
      float x = 0.0f;
      float y = 0.0f;
      float z = 0.0f;	
   };


   struct Vec3
   {
      float x = 0.0f;
      float y = 0.0f;
      float z = 0.0f;
   };
   Vec3 Vector(const float& x, const float& y, const float& z);


   class BlueIMU
   {

   public:



      BlueIMU(int gyroSampleRate);

      ///// SET FUNCTIONS
      void setGyroSampleRate( int rate);
      void setMagDeclination( float decDegrees);
      void setHeadingWithMag(const float& forward,const float& right, const float& down);
      void enableVelAndPosEstimation(int accelSampleRate);

      ///// INPUT FUNCTIONS
      void inputGyro( const float& rollRight, const float& pitchUp, const float& yawRight);
      void inputAccel(const float& forward,   const float& right,   const float& down);
      void inputMag(  const float& forward,   const float& right,   const float& down);
      void velocityCorrection(const float& forward, const float& right, const float& down, const float& weight);
      void positionCorrection(const float& north,   const float& east,  const float& down, const float& weight);

      ///// OPERATION FUNCTIONS
      void computeDriftCompensation();


      ///// OUTPUT FUNCTIONS
      float rollRadians();
      float pitchRadians();
      float yawRadians();

      //float rollSine();
      //float pitchSine();

      ///// 



      //private:



      ///// VARIABLES
      Quat  attitudeEstimate,
            rotationDelta;

      Vec3  velocity, position;

      Vec3  accelAccumulator,
            magAccumulator,
            correctionVectorBody;

      Vec3  magneticNorth = Vector(1.0f, 0.0f, 0.0f);

      unsigned long  GyroSamplePeriod,
                     AccelSamplePeriod,
                     MagSamplePeriod;

      float magDeclination,
            accelDeltaTime;

      bool  newDriftCorrection,
            newAccelData,
            newMagData;

      bool integrateVelAndPos = false;
      
      uint8_t normalizeQuat = 0;

   };
	

   /*
    *	Blue IMU: MATH FUNCTIONS
    */

   ///// CONSTANTS
   const Vec3 VERTICAL = Vector(0.0f, 0.0f, -1.0f);
   const Vec3 ZEROES   = Vector(0.0f, 0.0f, 0.0f);
   const float Deg2Rad = 0.017532925f;
   const float Rad2Deg = 57.29577951f;
          
   inline Quat Mul(const Quat& lhs, const Quat& rhs)  // multiply: Quat * Quat --16 mult, 12 add/sub  (266us)
   {
      Quat a;
      a.x = lhs.w * rhs.x + lhs.z * rhs.y - lhs.y * rhs.z + lhs.x * rhs.w;  
      a.y = lhs.w * rhs.y + lhs.x * rhs.z + lhs.y * rhs.w - lhs.z * rhs.x;
      a.z = lhs.y * rhs.x - lhs.x * rhs.y + lhs.w * rhs.z + lhs.z * rhs.w;
      a.w = lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z;
      return a;
   }


   inline Quat Mul(const Quat& lhs, const Vec3& rhs)  // multiply: Quat * [0,vec3] --12mult, 3 add, 6 sub
   {
      Quat a;
      a.x =  lhs.w * rhs.x + lhs.z * rhs.y - lhs.y * rhs.z;  
      a.y =  lhs.w * rhs.y + lhs.x * rhs.z - lhs.z * rhs.x;
      a.z =  lhs.y * rhs.x - lhs.x * rhs.y + lhs.w * rhs.z;
      a.w = -lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z;
      return a;
   }

   inline Quat Mul(Quat a, const float& b) // multiply: quat * float
   {
      a.w *= b;
      a.x *= b;
      a.y *= b;
      a.z *= b;
      return a;
   }

   inline Vec3 Mul(Vec3 a, const float& b) // multiply: vec3 * float
   {
      a.x *= b;
      a.y *= b;
      a.z *= b;
      return a;
   }


   inline Vec3 Mul(const float& b, const Vec3& a) // multiply: float * vec3
   {
      return Mul(a, b);
   }


   inline Vec3 Sum(Vec3 a, const Vec3& b)
   {
      a.x += b.x;
      a.y += b.y;
      a.z += b.z;
      return a;
   }


   inline float Magnitude(const Quat& a)
   {
      return sqrt(a.w*a.w + a.x*a.x + a.y*a.y + a.z*a.z);
   }


   inline float Magnitude(const Vec3& a)
   {
      return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
   }


   inline Vec3 CrossProd(const Vec3& L, const Vec3& R) // cross product of 3D vectors 
   { 
      Vec3 a;
      a.x = L.y * R.z - L.z * R.y;
      a.y = L.z * R.x - L.x * R.z;
      a.z = L.x * R.y - L.y * R.x;
      return a;	// 76us
   }


   inline Quat Quaternion(const Vec3& w, const unsigned long& t)   // (angular vel vector[rad/s], time interval[us]) -- Small angle approximation
   {  
      float dT_2;
      Quat  a;    

      dT_2 = float(t) * 0.0000005f; // time in seconds & divided in half for theta/2 computations
      a.x  = w.x * dT_2;
      a.y  = w.y * dT_2;
      a.z  = w.z * dT_2;
      a.w  = 1.0f - 0.5f * (a.x * a.x + a.y * a.y + a.z * a.z);
      return a;   // time = 116us + mult = 362us      (REF: RotMatrix = 588us)
   }


   inline Quat Quaternion(const Vec3& w)    // (angle vector[rad])    --Large Rotation Quaternion
   {  
      Quat a;

      float vMag = Magnitude(w);

      if(vMag < 1.0e-6) return a;  // terminate early if magnitude is negligible (a is zero rotation quat 0,0,0,1)

      float theta_2 = vMag * 0.5f;            // rotation angle divided by 2
      float Sin_Mag = sin(theta_2) / vMag;    // computation minimization

      a.x = w.x * Sin_Mag;
      a.y = w.y * Sin_Mag;
      a.z = w.z * Sin_Mag;
      a.w = cos(theta_2);

      return a;        // time = 390us + mult  = 636
   }


   inline Vec3 Vector(const float& x, const float& y, const float& z) // 3x float to vector
   { 
      Vec3 a;
      a.x = x;
      a.y = y;
      a.z = z;
      return a;
   }


   inline Vec3 Rotate(const Quat& q, const Vec3& v)	// Vector rotated by a Quaternion (matches V^ = Matrix * V)
   {	
      // v + 2*r X (r X v + q.w*v) -- https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
      Vec3 r; // vector r is the three imaginary coefficients of quaternion q
      r.x = q.x;
      r.y = q.y;
      r.z = q.z;
      return Sum(v, CrossProd(Sum(r, r), Sum(CrossProd(r, v), Mul(q.w, v))));	// 296us	
   }


   inline Vec3 Rotate(const Vec3& v, const Quat& q)	// Vector rotated by a Quaternion (matches V^ = V * Matrix)
   {	
      // v + 2*r X (r X v + q.w*v) -- https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
      Vec3 r;		// vector r is the three imaginary coefficients of quaternion q
      r.x = -q.x;	// reverse signs to change direction of rotation
      r.y = -q.y;
      r.z = -q.z;
      return Sum(v, CrossProd(Sum(r, r), Sum(CrossProd(r, v), Mul(q.w, v)))); // 296us 
   }


   /*
    *	Blue UMU - Data Display Functions
    *	
    *	SERIAL_PORT should be defined in the project
    *	example: #define SERIAL_PORT Serial1
    */

   /* 
   void display(const Vec3& v)
   {
      #ifdef SERIAL_PORT
         String outputBuffer;

         outputBuffer  = "X: ";
         outputBuffer += String(v.x, 4);
         outputBuffer += "   Y: ";
         outputBuffer += String(v.y, 4);
         outputBuffer += "   Z: ";
         outputBuffer += String(v.z, 4);
         outputBuffer += '\n';

         SERIAL_PORT.print(outputBuffer);
      #endif     
   }


   void display(const Quat& q)
   {
      #ifdef SERIAL_PORT
         String outputBuffer;

         outputBuffer  = "X: ";
         outputBuffer += String(q.x, 4);
         outputBuffer += "  Y: ";
         outputBuffer += String(q.y, 4);
         outputBuffer += "  Z: ";
         outputBuffer += String(q.z, 4);
         outputBuffer += "  W: ";
         outputBuffer += String(q.w, 4);
         outputBuffer += '\n';

         SERIAL_PORT.print(outputBuffer);
      #endif  
   } 
   */
#endif