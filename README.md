# float_to_decimal
A simple implementation of float to decimal in C that almost anyone can understand.

This was written less as a library, and more as a reference implenetation for anyone who wishes to copy-paste into their own codebase. However I did make some changea to make it esier to use as a standalone library.

Add `#define FLOAT_TO_DECIMAL_IMPLEMENTATION` before `#include`ing this file in *one* C or C++ file to create the implementation.

The main functions are:
```C
size_t FloatToDecimal(uint64_t mantisa, int32_t exponent, char* buffer, int32_t maxPrecision);
size_t Float32ToDecimal(float value, char* buffer, int32_t maxPrecision);
size_t Float64ToDecimal(double value, char* buffer, int32_t maxPrecision);
```

The return value is the length of the string. It is the caller's responsiblity to allocate enough memory.\
`maxPrecision` is the max number of digits after the radix.\
Passing -1 to `maxPrecision` argument will result in all the digits I know how to produce.\
Correct rounding is not implenmented.\
No null-terminator is added at the end.
