# float_to_decimal
A simple implementation of float to decimal in C-style C++ that almost anyone can understand.

The main functions are:
```C
size_t FloatToDecimal(uint64_t mantisa, int32_t exponent, int32_t precision, char* buffer);
size_t Float32ToDecimal(float value, int32_t precision, char* buffer);
size_t Float64ToDecimal(double value, int32_t precision, char* buffer);
```

The return value is the length of the string. It is the caller's respobsiblity to allocate enough memory.\
Passing -1 to `precision` argument will result in all the digits I know how to produce.\
Correct rounding is not implenmented.\
No null-terminator is added at the end.
