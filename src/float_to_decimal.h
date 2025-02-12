#if defined(__GNUC__)
#  define COMPILER_GCC 1
#elif defined(__clang__)
#  define COMPILER_CLANG 1
#elif defined(_MSC_VER)
#  define COMPILER_MSVC 1
#endif

#if defined(COMPILER_MSVC)
#  include <tmmintrin.h>
#  include <wmmintrin.h>
#  include <intrin.h>
#endif

#include <memory.h>
#include <immintrin.h>
#include <stdint.h>

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
# define __debugbreak() __builtin_trap()
#endif

#if defined(DEBUG)
#  define ASSERT(expression) 	do{if(!(expression)) {__debugbreak();}}while(0)
#else
#  define ASSERT(expression)	(void)(expression)
#endif

static uint64_t log10(uint64_t n) {
	if (n < 10000000000) {
		// 1..10
		return n < 100000
			? n < 1000
				? n < 100
					? n < 10 ? 1 : 2
					: 3
			
				: n < 10000 ? 4 : 5

			: n < 100000000
				? n < 10000000
					? n < 1000000 ? 6 : 7
					: 8

				: n < 1000000000 ? 9 : 10;
	}
	else {
		// 11..20
		return n < 1000000000000000ull
			? n < 10000000000000ull
				? n < 1000000000000ull
					? n < 100000000000ull ? 11 : 12
					: 13
			
				: n < 100000000000000ull ? 14 : 15

			: n < 1000000000000000000ull
				? n < 100000000000000000ull
					? n < 10000000000000000ull ? 16 : 17
					: 18

				: n < 10000000000000000000ull ? 19 : 20;
	}
}

// NOTE: don't use if divisor is known at compile time
static inline uint64_t udiv(uint64_t high, uint64_t low, uint64_t divisor, uint64_t* remainder) {
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	uint64_t quotient;
	uint64_t rem;

	asm("divq %4"
		: "=a"(quotient), "=d"(rem)
		: "d"(high), "a"(low), "r"(divisor));

	*remainder = rem;

	return quotient;
#elif defined(COMPILER_MSVC)
	return _udiv128(high, low, divisor, remainder);
#endif
}

// NOTE: don't use if divisor is known at compile time
static inline uint64_t udiv(uint64_t high, uint64_t low, uint64_t divisor, uint64_t* remainder, uint64_t* high_quotient) {
	uint64_t high_r = high % divisor;
	*high_quotient = high / divisor;
	return udiv(high_r, low, divisor, remainder);
}

static uint64_t udiv5(uint64_t high, uint64_t low, uint64_t* remainder, uint64_t* high_quotient) {
    // Use 2^64 = 1 mod 5 to compute remainder:
    // (high * 2^64 + low) 
    // = high + low mod 5
    // = high + low + carry mod 5 (the high + low might overflow, carry = 2^64 = 1 mod 5)
    //
    // Write high * 2^64 + low - remainder = highSub * 2^64 + lowSub
    // Note that highSub + lowSub = highSub * 2^64 + lowSub = 0 mod 5
    // Compute low bits of the divide:
    // lowSub * (2^66+1)/5 (division is exact!)
    // = lowSub/5 * 2^66 + lowSub/5
    // = lowSub*4/5 * 2^64 + lowSub/5
    // = highSub/5 * 2^64 + lowSub/5 mod 2^64 (lowSub * 4 = -lowSub = highSub mod 5)
    // = (highSub * 2^64 + lowSub)/5 mod 2^64
    // = flowor((high * 2^64 + low)/5) mod 2^64

    uint64_t merged;
    unsigned char carry = _addcarry_u64(0, high, low, &merged);
    _addcarry_u64(carry, merged, 0, &merged);
    uint64_t rem = merged % 5;
    uint64_t lowSub = low - rem;
    uint64_t low_quotient = lowSub * 14757395258967641293ull;
    *high_quotient = high / 5;

    *remainder = rem;

    return low_quotient;
}

static uint64_t udiv10(uint64_t high, uint64_t low, uint64_t* remainder, uint64_t* quotient_high) {
    // low2 = __shighftright128(low, high, 1);
    uint64_t low2 = (high << 63) | (low >> 1);
    uint64_t high2 = (high >> 1);
   
    uint64_t merged;
    unsigned char carry = _addcarry_u64(0, high2, low2, &merged);
    _addcarry_u64(carry, merged, 0, &merged);

    uint64_t rem = merged % 5;
    uint64_t lowSub = low2 - rem;
    uint64_t low_quotient = lowSub * 14757395258967641293ull;

    *remainder = low - 2*lowSub;
    *quotient_high = high / 10;

    return low_quotient;
}

static uint64_t umul(uint64_t a, uint64_t b, uint64_t *high_c) {
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	__uint128_t product = (__uint128_t)a * (__uint128_t)b;
	*high_c = product >> 64;
	return product;
#elif defined(COMPILER_MSVC)
	return _umul128(a, b, high_c);
#endif
}

static uint64_t umul(uint64_t high_a, uint64_t low_a, uint64_t b, uint64_t *high_c) {
	uint64_t low_c = umul(low_a, b, high_c);
	*high_c += high_a*b;
	return low_c;
}

inline static int64_t HighBit(uint64_t value, int64_t onZero) {
	if (value == 0) return onZero;

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	return 63 - __builtin_clzll(value);
#elif defined(COMPILER_MSVC)
	unsigned long index;
	_BitScanReverse64(&index, value);
	return (int64_t)index;
#endif
}

inline static int64_t LowBit(uint64_t value) {
	if (value == 0) return 64;

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	return __builtin_ctzll(value);
#elif defined(COMPILER_MSVC)
	unsigned long index;
	_BitScanForward64(&index, value);
	return (int64_t)index;
#endif
}

size_t UnsignedToDecimal(uint64_t number, char* str) {
	size_t numberOfDigits = log10(number);
	size_t index = numberOfDigits;

	do {
		index--;
		str[index] = '0' + (number % 10);
		number /= 10;
	} while (index);

	return numberOfDigits;
}

// TODO: I don't love this implementation
size_t UnsignedToDecimal(uint64_t high, uint64_t low, char* str) {
	if (high == 0)
		return UnsignedToDecimal(low, str);

	// get the last decimal digit
	uint64_t last;
	low = udiv10(high, low, &last, &high);

	char* ptr = str;
	if (high == 0) {
		ptr += UnsignedToDecimal(low, ptr);
	}
	else {
		// TODO: division by a constant should not be done with div
		uint64_t remainder;
		uint64_t quotient = udiv(high, low, 10000000000000000000ull, &remainder);
		ptr += UnsignedToDecimal(quotient, ptr);
		// TODO: I might want to add a pad parameter to UnsignedToDecimal
		size_t length = UnsignedToDecimal(remainder, ptr);
		if (length < 19) {
			size_t diff = 19 - length;
			memmove(ptr + diff, ptr, length);
			memset(ptr, '0', diff);
		}
		ptr += 19;
	}
	*(ptr++) = (char)last + '0';
	return ptr - str;
}

static int32_t bitIndexToExp[53] = {
	55, 54, 53, 53, 52, 52, 52, 51, 51, 50, 50, 
	49, 49, 49, 48, 48, 47, 47, 46, 46, 46, 45, 45, 
	44, 44, 43, 43, 43, 42, 42, 41, 41, 40, 40, 40, 
	39, 39, 38, 38, 37, 37, 37, 36, 36, 35, 35, 
	34, 34, 34, 33, 33, 32, 32
};

static struct {uint64_t high, low;} powersOf5[24] = {
	{0x4ee, 0x2d6d415b85acef81},
	{0x18a6, 0xe32246c99c60ad85},
	{0x7b42, 0x6fab61f00de36399},
	{0x2684c, 0x2e58e9b04570f1fd},
	{0xc097c, 0xe7bc90715b34b9f1},
	{0x3c2f70, 0x86aed236c807a1b5},
	{0x12ced32, 0xa16a1b11e8262889},
	{0x5e0a1fd, 0x2712875988becaad},
	{0x1d6329f1, 0xc35ca4bfabb9f561},
	{0x92efd1b8, 0xd0cf37be5aa1cae5},
	{0x2deaf189c, 0x140c16b7c528f679},
	{0xe596b7b0c, 0x643c7196d9ccd05d},
	{0x47bf19673d, 0xf52e37f2410011d1},
	{0x166bb7f0435, 0xc9e717bb45005915},
	{0x701a97b150c, 0xf18376a85901bd69},
	{0x23084f676940, 0xb7915149bd08b30d},
	{0xaf298d050e43, 0x95d69670b12b7f41},
	{0x36bcfc1194751, 0xed30f03375d97c45},
	{0x111b0ec57e6499, 0xa1f4b1014d3f6d59},
	{0x558749db77f700, 0x29c77506823d22bd},
	{0x1aba4714957d300, 0xd0e549208b31adb1},
	{0x85a36366eb71f04, 0x147a6da2b7f86475},
	{0x29c30f1029939b14, 0x6664242d97d9f649},
	{0xd0cf4b50cfe20765, 0xfff4b4e3f741cf6d},
};

size_t FloatToDecimal(uint64_t m2, int32_t e2, int32_t precision, char* buffer) {
	if (m2 == 0) {
		buffer[0] = '0';
		buffer[1] = '.';
		buffer[2] = '0';
		return 3;
	}

	uint32_t p = (uint32_t)precision;
	char* ptr = buffer;

	int32_t trailingZeroes = (int32_t)LowBit(m2);
	m2 >>= trailingZeroes;
	e2 += trailingZeroes;

	int32_t highBit = (int32_t)HighBit(m2, 0);
	ASSERT(highBit <= 52);
	bool hasWhole = 0 <= e2 || -e2 <= highBit;
	if (hasWhole) {
		bool wholeFitsIn128 = highBit + e2 < 128;
		if (wholeFitsIn128) {
			uint64_t high = 0;
			uint64_t low = 0;
			if (e2 < 0) {
				low = m2 >> -e2;
			}
			else if (e2 == 0) {
				low = m2;
			}
			else /*(e2 > 0)*/ {
				low = (e2 < 64) ? m2 << e2 : 0;
				high = (e2 < 64) ? m2 >> (64 - e2) : m2 << (e2 - 64);
			}
			ptr += UnsignedToDecimal(high, low, ptr);
		}
		else {
			uint64_t low = 0;
			uint64_t high = m2 << (63 - highBit);
			uint32_t e10 = 0;

			for (int32_t i = 0; i < e2 + highBit - 127; i++) {
				if ((high & 0x8000000000000000) == 0) {
					high = (high << 1) | ((low & 0x8000000000000000) >> 63);
					low <<= 1;
				}
				else {
					uint64_t remainder;
					low = udiv5(high, low, &remainder, &high);
					e10++;
				}
			}

			size_t length = UnsignedToDecimal(high, low, ptr);
			uint32_t exp = e10 + (int32_t)length - 1;
			if (length - 2 > p) length = p + 2;
			memmove(ptr + 1, ptr, length);
			ptr[1] = '.';
			ptr += length + 1;
			while (*(ptr - 1) == '0') ptr--;

			*(ptr++) = 'e';
			ptr += UnsignedToDecimal(exp, ptr);
			return ptr - buffer;
		}
	} else {
		*(ptr++) = '0';
	}

	*(ptr++) = '.';

	bool hasFraction = e2 < 0;
	if (hasFraction) {
		int32_t denominatorExp = -e2;
		if (denominatorExp < 64) {
			uint64_t denominator = 1ull << denominatorExp;
			uint64_t mask = denominator - 1;
			uint64_t numerator = m2 & mask;

			for (uint32_t digits = 0; numerator && digits < p; digits++) {
				uint64_t high;
				numerator = umul(numerator, 10, &high);
				// digit = udiv(high, numerator, denominator, &numerator)
				uint64_t digit = numerator >> denominatorExp | ((high & mask) << denominatorExp);
				numerator = numerator & mask;
				*(ptr++) = (char)digit + '0';
			}

			// remove trailing zeroes
			while (*(ptr - 1) == '0') ptr--;
		}
		else {
			ASSERT(!hasWhole);
			ptr -= 2;

			// we need to either multiply m2 by 5, or shift it by 1, -e2 times
			// start by multiplying m2 by 5^n, where n can be determined by the high-bit
			int32_t e10 = bitIndexToExp[highBit];
			if (e10 > -e2) e10 = -e2;
			uint64_t high;
			uint64_t low = umul(powersOf5[e10 - 32].high, powersOf5[e10 - 32].low, m2, &high);
			
			// we already multipled by 5 e10 times, so we don't start with 0
			for (int32_t i = e10; i < -e2; i++) {
				// if we can safely multiply by 5, then this is our preference
				if (high < 3689348814741910323) {
					low = umul(high, low, 5, &high);
					e10++;
				}
				else {
					low >>= 1;
					low |= ((high & 1) << 63);
					high >>= 1;
				}
			}

			size_t length = UnsignedToDecimal(high, low, ptr);
			int32_t exp = e10 - (int32_t)length + 1;
			if (length - 2 > p) length = p + 2;
			memmove(ptr + 1, ptr, length);
			ptr[1] = '.';
			ptr += length + 1;
			while (*(ptr - 1) == '0') ptr--;

			*(ptr++) = 'e';
			*(ptr++) = '-';
			ptr += UnsignedToDecimal(exp, ptr);
		}
	}
	else {
		*(ptr++) = '0';
	}

	return ptr - buffer;
}

size_t Float32ToDecimal(uint32_t value, int32_t precision, char* buffer) {
	uint32_t exp  = (value & 0x7F800000) >> 23;
	uint32_t sign = (value & 0x80000000) >> 31;
	uint64_t significand = (value & 0x7FFFFF);

	if (exp == 0xFF) {
		if (significand) {
			buffer[0] = 'N';
			buffer[1] = 'a';
			buffer[2] = 'N';
			return 3;
		}
		
		if (sign) {
			buffer[0] = '-';
			buffer[1] = 'i';
			buffer[2] = 'n';
			buffer[3] = 'f';
			return 4;
		}

		else {
			buffer[0] = 'i';
			buffer[1] = 'n';
			buffer[2] = 'f';
			return 3;
		}
	}

	size_t length = 0;
	if (sign) {
		buffer[0] = '-';
		buffer++;
		length++;
	}

	if (exp == 0) {
		length += FloatToDecimal(significand, -149, precision, buffer);
	}
	else {
		length += FloatToDecimal(significand | 0x800000, (int32_t)exp - 127 - 23, precision, buffer);
	}

	return length;
}

size_t Float32ToDecimal(float value, int32_t precision, char* buffer) {
	union {float f; uint32_t u;} cvt;
	cvt.f = value;
	return Float32ToDecimal(cvt.u, precision, buffer);
}

size_t Float64ToDecimal(uint64_t value, int32_t precision, char* buffer) {
	uint32_t sign =        (value & 0x8000000000000000) >> 63;
	uint32_t exp  =        (value & 0x7FF0000000000000) >> 52;
	uint64_t significand = (value & 0x000FFFFFFFFFFFFF) >> 00;

	if (exp == 0x7FF) {
		if (significand) {
			buffer[0] = 'N';
			buffer[1] = 'a';
			buffer[2] = 'N';
			return 3;
		}
		
		if (sign) {
			buffer[0] = '-';
			buffer[1] = 'i';
			buffer[2] = 'n';
			buffer[3] = 'f';
			return 4;
		}

		else {
			buffer[0] = 'i';
			buffer[1] = 'n';
			buffer[2] = 'f';
			return 3;
		}
	}

	size_t length = 0;
	if (sign) {
		buffer[0] = '-';
		buffer++;
		length++;
	}

	if (exp == 0) {
		length += FloatToDecimal(significand, -1075, precision, buffer);
	}
	else {
		length += FloatToDecimal(significand | 0x10000000000000, (int32_t)exp - 1023 - 52, precision, buffer);
	}

	return length;
}

size_t Float64ToDecimal(double value, int32_t precision, char* buffer) {
	union {double f; uint64_t u;} cvt;
	cvt.f = value;
	return Float64ToDecimal(cvt.u, precision, buffer);
}