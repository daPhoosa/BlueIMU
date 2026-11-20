# BlueIMU

A lightweight, quaternion-based attitude estimator (AHRS) for Arduino. It integrates gyroscope rates into an attitude quaternion. It corrects the resulting drift using accelerometer readings for roll and pitch, and magnetometer readings for yaw. It can also integrate acceleration into a velocity and position estimate in the world frame, which you can correct with external fixes such as GPS.

All math is single-precision `float`, so it runs efficiently on MCUs with a single-precision FPU (Teensy 3.x/4.x, ESP32, STM32F4, and others). It also builds on AVR boards.

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Frames, units and sign conventions](#frames-units-and-sign-conventions)
- [API reference](#api-reference)
- [Euler angles and gimbal lock](#euler-angles-and-gimbal-lock)
- [Fast trig option](#fast-trig-option)
- [Limitations](#limitations)
- [Planned](#planned)
- [License](#license)

## Installation

Copy the `BlueIMU` folder (containing `BlueIMU.h` and `BlueIMU.cpp`) into your Arduino `libraries` directory, then:

```cpp
#include <BlueIMU.h>
```

For PlatformIO, place the folder under `lib/`.

## Quick start

```cpp
#include <BlueIMU.h>

const int GYRO_RATE_HZ  = 500;   // must match how often inputGyro() is actually called
const int ACCEL_RATE_HZ = 500;   // only needed for velocity/position estimation

BlueIMU imu(GYRO_RATE_HZ);

void setup()
{
   Serial.begin(115200);
   // ... initialize your sensors ...

   imu.setMagDeclination(-3.5f);              // degrees, east positive (look up your location)
   imu.enableVelAndPosEstimation(ACCEL_RATE_HZ);  // optional

   // Optional: snap the initial heading to the magnetometer instead of waiting for drift
   // correction to converge. Works best with the device roughly level and still.
   float mf, mr, md;
   readMag(mf, mr, md);                        // your sensor read
   imu.setHeadingWithMag(mf, mr, md);
}

void loop()
{
   static unsigned long lastGyro = micros();
   static uint8_t       count    = 0;

   if (micros() - lastGyro >= 1000000UL / GYRO_RATE_HZ)
   {
      lastGyro += 1000000UL / GYRO_RATE_HZ;

      float gx, gy, gz;  readGyro(gx, gy, gz);   // rad/s
      float ax, ay, az;  readAccel(ax, ay, az);  // G
      imu.inputGyro(gx, gy, gz);
      imu.inputAccel(ax, ay, az);

      if (++count >= 10)                          // e.g. 50 Hz drift correction
      {
         count = 0;
         float mf, mr, md;  readMag(mf, mr, md);
         imu.inputMag(mf, mr, md);
         imu.computeDriftCompensation();
      }
   }

   // Optional: blend in a GPS fix when one arrives (world frame, m/s and m)
   // imu.velocityCorrection(velN, velE, velD, 0.1f);
   // imu.positionCorrection(posN, posE, posD, 0.1f);

   float roll, pitch, yaw;
   imu.getEulerRadians(roll, pitch, yaw);
}
```

`readGyro`, `readAccel` and `readMag` stand in for your sensor driver. Before passing readings in, remap the sensor axes to the body frame described below.

## Frames, units and sign conventions

**Body frame** (sensor inputs): X = forward, Y = right, Z = down.

**World frame** (velocity, position, corrections): NED, meaning X = north, Y = east, Z = down.

| Input | Units | Notes |
|---|---|---|
| Gyro | rad/s | `rollRight`, `pitchUp`, `yawRight` are positive for right-wing-down, nose-up, and nose-right (clockwise seen from above). |
| Accel | G | Specific force, as accelerometers report it. Level and at rest, the reading is `(0, 0, -1)`. |
| Mag | any | Only the direction is used, so the units don't matter. |
| Declination | degrees | East positive. |
| Velocity / position | m/s, m | World frame (NED). |
| Euler outputs | radians | Aerospace ZYX order. Roll: right wing down positive. Pitch: nose up positive. Yaw: clockwise from true north positive. |

## API reference

### Setup

| Function | Description |
|---|---|
| `BlueIMU(int gyroSampleRate)` | Creates the estimator. `gyroSampleRate` is in Hz and sets the fixed time step used to integrate each gyro sample. |
| `setGyroSampleRate(int rate)` | Changes the gyro rate (Hz). |
| `setMagDeclination(float degrees)` | Sets the angle between magnetic and true north (east positive). The yaw outputs are then relative to true north. |
| `setHeadingWithMag(forward, right, down)` | Rotates the current estimate about the vertical so that its heading matches the given magnetometer reading, with declination applied. Roll and pitch are kept. Useful at startup. |
| `enableVelAndPosEstimation(int accelSampleRate)` | Turns on velocity and position integration. Each `inputAccel()` call then advances them by `1 / accelSampleRate` seconds. |

### Inputs

| Function | Description |
|---|---|
| `inputGyro(rollRight, pitchUp, yawRight)` | Integrates one gyro sample into the attitude. Call this at exactly the configured rate. If drift compensation was computed since the last call, it is applied in this step. |
| `inputAccel(forward, right, down)` | Accumulates the reading for the next drift correction. If enabled, it also integrates velocity and position, with gravity removed. |
| `inputMag(forward, right, down)` | Accumulates the horizontal component of the reading for the next yaw correction. |
| `velocityCorrection(north, east, down, weight)` | Blends an external velocity (m/s, NED) into the estimate: `v = (1 - weight)·v + weight·measured`, where `weight` is between 0 and 1. |
| `positionCorrection(north, east, down, weight)` | The same blend, for position (m, NED). |

### Drift correction

| Function | Description |
|---|---|
| `computeDriftCompensation()` | Turns the accel and mag readings accumulated since the last call into a correction rate, then clears the accumulators. The next `inputGyro()` applies it. Call it periodically, typically at 10–100 Hz. |

### Outputs

| Function | Description |
|---|---|
| `float rollRadians() const` | Roll angle. |
| `float pitchRadians() const` | Pitch angle. |
| `float yawRadians() const` | Yaw (heading) angle. |
| `void getEulerRadians(float& roll, float& pitch, float& yaw) const` | All three at once. This shares intermediate terms, so it's cheaper than calling the three functions individually. |

Multiply by `Rad2Deg` to convert to degrees.

The raw state is also accessible as public members: `attitudeEstimate` (a `Quat`, Hamilton convention, rotates body to world), `velocity` and `position` (`Vec3`, NED).

### Math helpers

`BlueIMU.h` also provides small `Quat` and `Vec3` types with inline helpers:

- `Vector(x, y, z)` builds a `Vec3`.
- `Sum`, `Mul`, `CrossProd` and `Magnitude` do vector and quaternion arithmetic.
- `Rotate(q, v)` rotates `v` from body to world. `Rotate(v, q)` rotates from world to body.
- `Quaternion(angleVec)` builds a rotation quaternion from a rotation vector.

Constants: `Deg2Rad`, `Rad2Deg`, `GRAVITY` (9.80665 m/s²), `VERTICAL`, `ZEROES`.

## Euler angles and gimbal lock

The Euler outputs are computed directly from the quaternion. They are designed to stay accurate and well behaved at every attitude:

- **Pitch** uses `atan2(sin, cos)` rather than `asin`. This keeps it accurate all the way to ±90° without any clamping.
- **No normalization required.** Every formula is written so that a quaternion that has drifted from unit length still gives exact angles.
- **Gimbal lock.** At pitch ±90°, roll and yaw describe the same rotation and can't be separated. Within about 0.057° of ±90° pitch, roll is reported as 0 and the combined rotation is reported in yaw. The output is always finite and deterministic, never NaN. The threshold is `GIMBAL_LOCK_COS` in `BlueIMU.cpp`.

## Fast trig option

Defining `BLUEIMU_FAST_TRIG` replaces `atan2f` with a polynomial approximation that has a maximum error of about 1e-5 rad (0.0006°). This speeds up the Euler outputs on MCUs where `atan2f` is slow.

Arduino compiles libraries separately from the sketch, so a `#define` in your sketch has no effect. Set it as a compiler flag instead. For example, in PlatformIO:

```ini
build_flags = -DBLUEIMU_FAST_TRIG
```

## Limitations

- **Fixed time steps.** Gyro integration assumes exactly `1 / gyroSampleRate` seconds between `inputGyro()` calls, and velocity integration assumes `1 / accelSampleRate` between `inputAccel()` calls. Timing jitter shows up directly as error.
- **No gain tuning or gyro-bias estimation.** The strength of the drift correction depends on how many accel and mag samples accumulate between `computeDriftCompensation()` calls. There is no integral term, so a constant gyro bias produces a small steady attitude offset.
- **Accelerometer assumed to measure gravity only.** Sustained linear acceleration, such as turns or speeding up and slowing down, pulls the roll and pitch estimate toward a false vertical.
- **Yaw correction uses only the horizontal component of the field.** Accuracy degrades at large tilt angles and with uncalibrated magnetometers. Hard-iron and soft-iron calibration is up to the caller.
- **Velocity and position drift quickly.** Without external corrections they are only useful over short intervals.

## Testing

Host-side unit tests live in `test/`. They use [doctest](https://github.com/doctest/doctest) (vendored as a single header) and a stub `Arduino.h`, so no board or Arduino install is needed:

```sh
make -C test          # builds and runs twice: libm atan2f, and -DBLUEIMU_FAST_TRIG
make -C test exact    # one variant only
make -C test ARGS="-ts='*gimbal*'"   # pass doctest options, e.g. filter by suite
```

The library is compiled with `-Wdouble-promotion -Wfloat-conversion -Werror`, so any accidental double-precision math fails the build.

## Planned

- Barometer input for vertical position.
- Native GPS input, as an alternative to calling `velocityCorrection()` and `positionCorrection()` manually.

## License

Copyright (c) 2016-2026 Phillip J Schmidt. Released under the [MIT License](LICENSE).

The vendored test framework in `test/vendor/` ([doctest](https://github.com/doctest/doctest)) is MIT licensed by its authors. See [test/vendor/doctest-LICENSE.txt](test/vendor/doctest-LICENSE.txt).
