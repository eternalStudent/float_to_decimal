#if defined(__GNUC__)
#  define COMPILER_GCC 1
#elif defined(__clang__)
#  define COMPILER_CLANG 1
#elif defined(_MSC_VER)
#  define COMPILER_MSVC 1
#endif

#if defined(COMPILER_MSVC)
#  include <intrin.h>
#endif

#include <string.h>
#include <stdint.h>

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
# define __debugbreak() __builtin_trap()
#endif

#if defined(DEBUG)
#  define ASSERT(expression) 	do{if(!(expression)) {__debugbreak();}}while(0)
#else
#  define ASSERT(expression)	(void)(expression)
#endif

static inline int64_t bits_findFirstSetBit(uint64_t value, int64_t onZero) {
	if (value == 0) return onZero;

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	return 63 - __builtin_clzll(value);
#elif defined(COMPILER_MSVC)
	unsigned long index;
	_BitScanReverse64(&index, value);
	return (int64_t)index;
#endif
}

static inline int64_t bits_findLastSetBit(uint64_t value) {
	if (value == 0) return 64;

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	return __builtin_ctzll(value);
#elif defined(COMPILER_MSVC)
	unsigned long index;
	_BitScanForward64(&index, value);
	return (int64_t)index;
#endif
}

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
typedef unsigned long long carry_t;
#elif defined(COMPILER_MSVC)
typedef unsigned char carry_t;
#endif

uint64_t uint64_addWithCarry(carry_t carry, uint64_t a, uint64_t b, carry_t* high_c) {
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	return __builtin_addcll(a, b, carry, high_c);
#elif defined(COMPILER_MSVC)
	uint64_t c;
	*high_c = _addcarry_u64(carry, a, b, &c);
	return c;
#endif
}

static inline uint64_t uint128_shiftright(uint64_t high, uint64_t low, unsigned char shift, uint64_t* high_shifted) {
	*high_shifted = high >> shift;
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	return (high << (64 - shift)) | (low >> shift);
#elif defined(COMPILER_MSVC)
	return __shiftright128(low, high, shift);
#endif
}

static inline uint64_t uint128_shiftleft(uint64_t high, uint64_t low, unsigned char shift, uint64_t* high_shifted) {
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	*high_shifted = (high << shift) | (low >> (64 - shift));
#elif defined(COMPILER_MSVC)
	*high_shifted =  __shiftleft128(low, high, shift);
#endif
	return low << shift;
}

static inline uint64_t uint128_div(uint64_t high, uint64_t low, uint64_t divisor, uint64_t* remainder) {
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	uint64_t quotient;
	uint64_t rem;

	asm("divq %4"
		: "=a"(quotient), "=d"(rem)
		: "d"(high), "a"(low), "r"(divisor)
		: "cc");

	*remainder = rem;

	return quotient;
#elif defined(COMPILER_MSVC)
	return _udiv128(high, low, divisor, remainder);
#endif
}

static inline uint64_t uint128_div(uint64_t high, uint64_t low, uint64_t divisor, uint64_t* remainder, uint64_t* high_quotient) {
	uint64_t high_r = high % divisor;
	*high_quotient = high / divisor;
	return uint128_div(high_r, low, divisor, remainder);
}

static inline uint64_t uint128_div5(uint64_t high, uint64_t low, uint64_t* remainder, uint64_t* high_quotient) {
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

	carry_t carry;
	uint64_t merged = uint64_addWithCarry(0, high, low, &carry);
	merged = uint64_addWithCarry(carry, merged, 0, &carry);
	uint64_t rem = merged % 5;
	uint64_t lowSub = low - rem;
	uint64_t low_quotient = lowSub * 14757395258967641293ull;
	*high_quotient = high / 5;

	*remainder = rem;

	return low_quotient;
}

static inline uint64_t uint128_div10(uint64_t high, uint64_t low, uint64_t* remainder, uint64_t* quotient_high) {
	uint64_t high2;
	uint64_t low2 = uint128_shiftright(high, low, 1, &high2);

	carry_t carry;
	uint64_t merged = uint64_addWithCarry(0, high2, low2, &carry);
	merged = uint64_addWithCarry(carry, merged, 0, &carry);

	uint64_t rem = merged % 5;
	uint64_t lowSub = low2 - rem;
	uint64_t low_quotient = lowSub * 14757395258967641293ull;

	*remainder = low - 2*lowSub;
	*quotient_high = high / 10;

	return low_quotient;
}

static inline uint64_t uint128_mul(uint64_t high_a, uint64_t low_a, uint64_t b, uint64_t *high_c) {
#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
	__uint128_t a = ((__uint128_t)high_a << 64) | (__uint128_t)low_a;
	__uint128_t c = (__uint128_t)a * (__uint128_t)b;
	*high_c = c >> 64;
	return c;
#elif defined(COMPILER_MSVC)
	uint64_t low_c = _umul128(low_a, b, high_c);
	*high_c += high_a*b;
	return low_c;
#endif
}

static inline uint64_t getNumberOfDecimalDigits(uint64_t n) {
	static uint64_t table[] = {
		9,
		99,
		999,
		9999,
		99999,
		999999,
		9999999,
		99999999,
		999999999,
		9999999999,
		99999999999,
		999999999999,
		9999999999999,
		99999999999999,
		999999999999999ull,
		9999999999999999ull,
		99999999999999999ull,
		999999999999999999ull,
		9999999999999999999ull
	};
	uint64_t y = (19 * bits_findFirstSetBit(n, 0) >> 6);
	y += n > table[y];
	return y + 1;
}

size_t UnsignedToDecimal(uint64_t number, char* str) {
	size_t numberOfDigits = getNumberOfDecimalDigits(number);
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

	/*
	 * 128bit unsigned integers can be up to 40 decimal digits.
	 * 64bit unsigned integer can be up to 20 decimal digits.
	 * You might think that you can just divide it into 2 parts:
	 * Print the most significant 20 decimal digits, and then the least significant 20 decimal digits.
	 * But the least significant 20 digits can actually add up to more than 64 bits.
	 * So I am splitting it up into 3 parts.
	 * This is not a great implementation, the uint128_div will produce an expansive div instruction
	 * Improvements suggestions are welcome.
	 */

	uint64_t last;
	low = uint128_div10(high, low, &last, &high);

	char* ptr = str;
	if (high == 0) {
		ptr += UnsignedToDecimal(low, ptr);
	}
	else {
		uint64_t remainder;
		uint64_t quotient = uint128_div(high, low, 10000000000000000000ull, &remainder);
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

size_t FloatToDecimal(uint64_t m2, int32_t e2, int32_t precision, char* buffer) {

	static int32_t exponentOf5[53] = {
		55, 54, 53, 53, 52, 52, 52, 51, 51, 50, 50, 
		49, 49, 49, 48, 48, 47, 47, 46, 46, 46, 45, 45, 
		44, 44, 43, 43, 43, 42, 42, 41, 41, 40, 40, 40, 
		39, 39, 38, 38, 37, 37, 37, 36, 36, 35, 35, 
		34, 34, 34, 33, 33, 32, 32
	};

	// starting with 5^32
	static struct {uint64_t high, low;} powerOf5[24] = {
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

	if (m2 == 0) {
		buffer[0] = '0';
		buffer[1] = '.';
		buffer[2] = '0';
		return 3;
	}

	uint32_t p = (uint32_t)precision;
	char* ptr = buffer;

	int32_t trailingZeroes = (int32_t)bits_findLastSetBit(m2);

	// simplify the fraction
	m2 >>= trailingZeroes;
	e2 += trailingZeroes;

	int32_t highBit = (int32_t)bits_findFirstSetBit(m2, 0); 
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
			/*
			 * We need to reduce the exponent of 2 to 0, 
			 * but we can only shift left `127 - highBit` times
			 * and still fit into 128 bits
			 * and we know that `highBit + e2 > 127`, so it's not enough.
			 *
			 * The trick here is that x*2^n = x*2^n * 5^(-m)*5^m
			 *                              = x*2^(n-m)*2^m * 5^(-m)*5^m
			 *                              = x*2^(n-m) * 5^(-m) * 2^m*5^m
			 *                              = x*2^(n-m)*5^(-m) * 10^m
			 *
			 * So each time we divide the mantisa by 5, 
			 * we reduce the exponent of 2 by 1
			 * and increase the exponent of 10 by 1
			 * but dividing by 5 reduces precision,
			 * so whenever possible we prefer to shift left.
			 */
			uint64_t low = 0;
			uint64_t high = m2 << (63 - highBit);
			uint32_t e10 = 0;

			for (int32_t i = 0; i < e2 + highBit - 127; i++) {
				if ((high & 0x8000000000000000) == 0) {
					low = uint128_shiftleft(high, low, 1, &high);
				}
				else {
					uint64_t remainder;
					low = uint128_div5(high, low, &remainder, &high);
					e10++;
				}
			}

			size_t length = UnsignedToDecimal(high, low, ptr);
			uint32_t exp = e10 + (uint32_t)length - 1;
			if (length - 2 > p) length = p + 2;
			memmove(ptr + 1, ptr, length);
			ptr[1] = '.';
			ptr += length + 1;

			// remove trailing zeroes
			while (*(ptr - 1) == '0') ptr--;

			*(ptr++) = 'e';
			ptr += UnsignedToDecimal(exp, ptr);

			// there's not going to be any fraction, so just return
			return ptr - buffer;
		}
	} else {
		*(ptr++) = '0';
	}

	*(ptr++) = '.';

	bool hasFraction = e2 < 0;
	if (hasFraction) {
		unsigned char denominatorExp = (unsigned char)-e2;
		if (denominatorExp < 64) {
			uint64_t denominator = 1ull << denominatorExp;
			uint64_t mask = denominator - 1;
			uint64_t numerator = m2 & mask; // remove the whole part

			// This is just long division
			for (uint32_t digits = 0; numerator && digits < p; digits++) {
				uint64_t high;
				numerator = uint128_mul(0, numerator, 10, &high);
				uint64_t digit = uint128_shiftright(high, numerator, denominatorExp, &high);
				numerator = numerator & mask;
				*(ptr++) = (char)digit + '0';
			}

			// remove trailing zeroes
			while (*(ptr - 1) == '0') ptr--;
		}
		else {
			ASSERT(!hasWhole);
			// oops we printed `0.` but we don't need it.
			ptr -= 2;

			/*
			 * We need to increase the exponent of 2 to 0.
			 *
			 * The trick here is that x*2^(-n) = x*(2^-n) * 5^m*5^(-m)
			 *                                 = x*2^(m-n)*2^m * 5^m*5^(-m)
			 *                                 = x*2^(m-n) * 5^m * 2^(-m)*5^(-m)
			 *                                 = x*2^(m-n)*5^m * 10^(-m)
			 *
			 * So each time we multiple the mantisa by 5, 
			 * we increase the exponent of 2 by 1
			 * and decrease the exponent of 10 by 1.
			 * We'll figure out how many times we can safely multiply by 5
			 * using the hight bit, and taking the worst case scenario.
			 * and then multiply by a 5 to the power of that many times.
			 * Whenever we can't multiply by 5 without overflowing, 
			 * we will shift right instead, which also increases the exponent of 2,
			 * but decreases precision.
			 */

			int32_t e10 = exponentOf5[highBit];
			if (e10 > -e2) e10 = -e2;
			uint64_t high;
			uint64_t low = uint128_mul(powerOf5[e10 - 32].high, powerOf5[e10 - 32].low, m2, &high);
			for (int32_t i = e10; i < -e2; i++) {
				if (high < 3689348814741910323ull) {
					low = uint128_mul(high, low, 5, &high);
					e10++;
				}
				else {
					low = uint128_shiftright(high, low, 1, &high);
				}
			}

			size_t length = UnsignedToDecimal(high, low, ptr);
			int32_t exp = e10 - (int32_t)length + 1;
			if (length - 2 > p) length = p + 2;
			memmove(ptr + 1, ptr, length);
			ptr[1] = '.';
			ptr += length + 1;

			// remove trailing zeroes
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