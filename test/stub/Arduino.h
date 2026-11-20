// Minimal host stand-in for the Arduino core: only what BlueIMU uses.
#pragma once
#include <math.h>
#include <stdint.h>

// Controllable clock: tests set/advance stubMicros() to drive BlueIMU's dt measurement.
inline unsigned long& stubMicros() { static unsigned long t = 0; return t; }
inline unsigned long  micros()     { return stubMicros(); }
