/*
   BlueIMU basic example

   Feeds gyro, accel and mag samples into BlueIMU and prints roll, pitch and yaw
   (degrees) ten times a second.

   The read*() functions below are placeholders that simulate a level, stationary
   sensor pointing north. Replace them with your sensor driver, remapping the
   sensor axes to BlueIMU's body frame (X forward, Y right, Z down).
 */

#include <BlueIMU.h>

BlueIMU imu;

// ---- Replace these with your sensor driver ----
void readGyro(float& rollRate, float& pitchRate, float& yawRate)  // rad/s
{
   rollRate = 0.0f;  pitchRate = 0.0f;  yawRate = 0.0f;
}

void readAccel(float& forward, float& right, float& down)          // G
{
   forward = 0.0f;  right = 0.0f;  down = -1.0f;                   // level and at rest
}

void readMag(float& forward, float& right, float& down)            // any units
{
   forward = 0.2f;  right = 0.0f;  down = 0.45f;                   // facing magnetic north
}
// -----------------------------------------------

void setup()
{
   Serial.begin(115200);

   imu.setMagDeclination(0.0f);       // degrees, east positive: look up your location
   imu.enableVelAndPosEstimation();   // optional

   float mf, mr, md;
   readMag(mf, mr, md);
   imu.setHeadingWithMag(mf, mr, md); // start with the correct heading
}

void loop()
{
   static unsigned long lastSample     = 0;
   static unsigned long lastCorrection = 0;
   static unsigned long lastPrint      = 0;
   unsigned long now = micros();

   // Sample at ~500 Hz. The rate doesn't need to be exact: BlueIMU measures the
   // interval between samples itself. With real sensors, read when data is ready.
   if (now - lastSample >= 2000UL)
   {
      lastSample = now;

      float gx, gy, gz;  readGyro(gx, gy, gz);
      float ax, ay, az;  readAccel(ax, ay, az);
      float mf, mr, md;  readMag(mf, mr, md);

      imu.inputGyro(gx, gy, gz);
      imu.inputAccel(ax, ay, az);
      imu.inputMag(mf, mr, md);
   }

   if (now - lastCorrection >= 20000UL)   // 50 Hz drift correction
   {
      lastCorrection = now;
      imu.computeDriftCompensation();
   }

   if (now - lastPrint >= 100000UL)       // 10 Hz output
   {
      lastPrint = now;

      float roll, pitch, yaw;
      imu.getEulerRadians(roll, pitch, yaw);

      Serial.print(F("roll "));     Serial.print(roll  * Rad2Deg, 2);
      Serial.print(F("  pitch "));  Serial.print(pitch * Rad2Deg, 2);
      Serial.print(F("  yaw "));    Serial.println(yaw * Rad2Deg, 2);
   }
}
