#include <iostream>
#include <array>
#include <gcem.hpp>
#include "pretty_types.h"

/* the fibBase data type will be used as the main and maximum allowed data type for
 * the slow fibonacci encoding algorithm, which is more useful/simpler/faster
 * for smaller numbers (like until 2^16)
 *
 * the fibBaseHalf data type is needed for some optimizations in the algorithms, so that we don't
 * have to figure out hacky workarounds just specify it here please!
 *
 * Sadly I couldn't figure out how to get the preprocessor to know the bit size of the type,
 * so you also need to specify the size in bits in FIB_BASE_BITS, this is needed for, again,
 * simplicity and optimizations
 *
 * making this larger, for example uint128, doesn't make sense, but I guess who will stop ya?
*/

typedef uint16 fibBase;
constexpr fibBase MAX_FIBBASE = std::numeric_limits<fibBase>::max();
constexpr uint32 fibBaseBits = sizeof(fibBase) * 8;

constexpr long double phi = (1.0L + gcem::sqrt(5.0L)) / 2.0L;
constexpr long double Phi = 1.0L - phi;

constexpr uint64 FastFibonacci(uint64 i) {
	return (gcem::pow(phi, i) - gcem::pow(Phi, i)) / gcem::sqrt(5.0L);
}

/* Function that returns the index of the largest fibonacci number that's less than n
 * (also the length of n encoded with fibonacci)
 *
 * Can be derived using Binet's representation and the fact that (1 - phi)^n -> 0 as n -> inf
 * And actually even from the first fibonacci numbers, it barely affects the value, and simply
 * Rounding we get the correct fibonacci numbers. So instead (phi^n + (1 - phi)^n)/sqrt(5) we have
 * Simply phi^n/sqrt(5) = ((1 + sqrt(5))/2)^n * 1/sqrt(5). Solving the inequality
 * ((1 + sqrt(5))/2)^n * 1/sqrt(5) <= N where N is the function argument, we get n - the index
 * of the largest fibonacci number that'll fit us i.e. be still smaller or equal to N.
 *
 * However because of this "rounding" we need to make a const table for small numbers, because
 * they won't be found correctly by this formula
 */
constexpr uint32 small_number_fitter[20] =
	{0, 2, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7};
constexpr uint32 MaxFittingFibonacci(uint64 n) {
		if (n < 20)
			return small_number_fitter[n];
		constexpr long double sqrt5 = gcem::sqrt(5.0L);
		constexpr long double denominator = gcem::log(1.0L + sqrt5) - gcem::log(2.0L);
		// add some epsilon to escape invalid answers for fibonacci numbers
		// higher than the ones stored in the small fitter
		return gcem::log(sqrt5 * static_cast<long double>(n)) / denominator + 0.001L;
}

/* Pretty much just a trunk for the last function, but handles 0 differently
 * because the encoding of 0 is nothing (actually just the '1' which is placed after every code)
 */
constexpr inline uint32 FibonacciEncodingLength(uint64 n) {
	if (n != 0) [[likely]] {
		return MaxFittingFibonacci(n) - 1;
	}
	return 0;
}

/* Build precomputed at compile time fibonacci table of size n
 *
 * We later use this in order to build the base encoding table for small values
 * As well as for simply converting small numbers to fibonacci encoding
 */
template <size_t n>
constexpr std::array<fibBase, n> BuildPrecomputedTable() {
	std::array<fibBase, n> result = {1};
	if (n > 1)
		result[1] = 2;
	for (uint32 i = 2; i < n; i++)
		result[i] = result[i - 1] + result[i - 2];
	return result;
}

// Select the "best" (actually not the best, of course, but a better one)
// fibonacci number for encoding numbers which are < 2^64 if we have fibBase available
//
constexpr uint32 SelectBestFibonacci() {
	// the maximum amount of bits we will need to store the max number in fibonacci encoding
	constexpr uint32 max_needed_fibonacci = MaxFittingFibonacci(std::numeric_limits<uint64>::max()) - 1;
	// index of the maximum fibonacci number which we can encode using fibBase
	constexpr uint32 max_available_fibonacci = MaxFittingFibonacci(MAX_FIBBASE) - 1;
	// the number of blocks we would need if each block would be max_av_fib bits long
	constexpr uint32 needed_blocks = max_needed_fibonacci / max_available_fibonacci;
	// the index of a fibonacci number which would still need needed_blocks amount of blocks,
	// but we would require less to store each table (for all the different shifts which are multiples of the block size)
	constexpr uint32 better_fit_index = (max_needed_fibonacci + needed_blocks - 1) / needed_blocks;
	return better_fit_index;
}

// Build a precomputed table with enough size to hold all fibonacci numbers that are less
// than the maximum number representable by our base type

// constexpr uint32 precalc_table_size = MaxFittingFibonacci(MAX_FIBBASE) - 1;

// Actually however we store the index of the maximum usable fibonacci number here, i.e. 16
// in the case of using an int16 as the base type. This is due to the fact that we encode
// the numbers in blocks of this size, and it is more optimal to store tables with sizes equal to
// some power of 2 than something like 25, which we would get if we stored numbers up to the maximal needed.
constexpr uint32 maxNeededFibonacci = MaxFittingFibonacci(MAX_FIBBASE) - 1;
constexpr uint32 fibBlockSize = fibBaseBits;
constexpr uint32 fibBlocksNeeded = maxNeededFibonacci / fibBlockSize;
constexpr fibBase precalc_table_limit = FastFibonacci(fibBlockSize + 2);
constexpr std::array<fibBase, fibBlockSize> precalc_fibonacci_table = BuildPrecomputedTable<fibBlockSize>();

// Length of the maximum code we can get, we only use up to uint64 so this will be enough for everything
// We will store the codes simply as bitsets, because there's no point in making some smart structure for them
// The overhead would just be too much, comparing to storing around 12 bytes
constexpr uint32 max_needed_fibonacci_index = MaxFittingFibonacci(std::numeric_limits<uint64>::max());

/* Lightweight class for storing a bitset of fixed size.
 *
 * This stores the bits in network order (if the system is little endian, otherwise we will get fucked...)
 * so that after we encode all the needed data we can just spew it out into the file instantly.
 *
 * It's easy to change this to a template class, but currently there's no need, cause we have a constant
 * maximum code length, and the overhead of a more complex class is bigger than just storing the extra bytes
 */
struct LightBitset {
	// the max length is the index of the maximum fibonacci number we can use, rounded up to 8-bit byte size
	static const uint32 LIGHT_BITSET_SIZE = ((max_needed_fibonacci_index + 7) >> 3);
	uint8 bits_[LIGHT_BITSET_SIZE];
	uint8 max_index_;

	constexpr LightBitset() {
		max_index_ = 0;
		Clear();
	}

	~LightBitset() = default;

	/* Basic functions for manipulating instances of this class */

	constexpr void Set(uint8 index) {
		bits_[(index >> 3)] |= (1 << (7 - (index & 0b111)));
		if (index + 1 > max_index_)
			max_index_ = index + 1;
	}

	constexpr bool Get(uint8 index) {
		return (bits_[(index >> 3)] & (1 << (7 - (index & 0b111)))) != 0;
	}

	constexpr void Clear() {
		for (uint32 i = 0; i < LIGHT_BITSET_SIZE; i++)
			bits_[i] = 0;
	}

	/* Function for merging this bitset with another one
	 *
	 * Pretty much just sets all the bits present in the other bitset
	 * Merges by block, not by bit. This is correct even if the last block in the other bitset
	 * isn't completely filled, because in that case the values will be 0 anyways
	 */
	LightBitset MergeWith(LightBitset &other) {
		LightBitset result;
		for (uint32 i = 0; i < ((other.max_index_ + 7) >> 3); i++) {
			result.bits_[i] = (bits_[i] | other.bits_[i]);
		}
		return result;
	}
};

/* Function for converting small numbers (up to MAX_FIBBASE) relatively easily
 *
 * This WILL NOT encode larger numbers correctly, because it uses the precomputed fibonacci number table
 * where all the numbers are lower than MAX_FIBBASE. This function should not be used to build large
 * consecutive tables, because you would have to encode every number using it, resulting in
 * O(n * log(precalc_table_size)) which is technically O(n) but much more inefficient
 */
constexpr LightBitset EncodeSmallNumber(fibBase n) {
	LightBitset result;
	int32 last = fibBlockSize;
	while (n > 0) {
		int left = 0, right = last;
		while (left < right) {
			int64 mid = (left + right) >> 1;
			if (precalc_fibonacci_table[mid] > n)
				right = mid;
			else
				left = mid + 1;
		}
		last = right;
		n -= precalc_fibonacci_table[right - 1];
		result.Set(right - 1);
	}
	return result;
}

/* This generates a code table for small numbers (lower than MAX_FIBBASE)
 *
 * Actually this can be calculated more efficiently if we take into account the fact
 * that we can use our previous results for encoding the next numbers:
 * 0 is ''
 * 1 is '1'
 * 2 is '01'
 * 3 is '001'
 * 4 is '101'
 * 5 is '0001'
 * 6 is '1001'
 * 7 is '0101'
 * 8 is '00001' and so on...
 * As you can see, we can encode all the numbers which have i'th fib. number as the maximum fitting fib number
 * by using all the previously encoded number UP UNTIL the previous fib. number (so for encoding numbers which
 * have 5 as max fitting fib. number, we use 1, 2, because 3 was previous fib number and we skip it, and we get
 * the encodings for 6 and 7) and then just append a '1' to the end, signifying this number. This can be easily
 * used to encode numbers in clear O(n) instead of amortized O(n) (which is due to using up lots of floating-point
 * arithmetic).
 * And since we even have a precalc'd table, we can do this really quite fast and get a quickly-compiling constexpr table.
 *
 * However it doesn't matter for now, so currently it just uses the simple encoding algorithm.
 */
template<size_t n>
constexpr std::array<LightBitset, n> GenerateCodeTable() {
	std::array<LightBitset, n> result = {};
	for (fibBase i = 0; i < n; i++) {
		result[i] = EncodeSmallNumber(i);
	}
	return result;
}

constexpr std::array<LightBitset, precalc_table_limit> precalc_code_table = GenerateCodeTable<precalc_table_limit>();



int main() {
	std::cout << max_needed_fibonacci_index << std::endl;
	std::cout << LightBitset::LIGHT_BITSET_SIZE << std::endl;
	std::cout << precalc_table_limit << std::endl;

	for (uint32 i = 0; i < fibBlockSize; i++) {
		std::cout << i << ": " << precalc_fibonacci_table[i] << std::endl;
	}
//	while (true) {
//		int i;
//		std::cin >> i;
//		std::cout << MaxFittingFibonacci(i) << " " << FibonacciEncodingLength(i) << std::endl;
//	}

	for (uint32 i = 0; i < 100; i++) {
		LightBitset encoded = precalc_code_table[i];
		std::cout << i << " ";
		for (uint32 j = 0; j < encoded.max_index_; j++) {
			std::cout << (encoded.Get(j) ? '1' : '0');
		}
		std::cout << '1' << std::endl;
	}

	return 0;
}
