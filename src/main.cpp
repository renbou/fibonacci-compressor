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

// Length of the maximum code we can get, we only use up to uint64 so this will be enough for everything
// We will store the codes simply as bitsets, because there's no point in making some smart structure for them
// The overhead would just be too much, comparing to storing around 12 bytes
constexpr uint32 max_needed_fibonacci_index = MaxFittingFibonacci(std::numeric_limits<uint64>::max());

// We store the index of the maximum usable fibonacci number here, i.e. 16
// in the case of using an int16 as the base type. This is due to the fact that we encode
// the numbers in blocks of this size, and it is more optimal to store tables with sizes equal to
// some power of 2 than something like 25, which we would get if we stored numbers up to the maximal needed.
constexpr uint32 fibonacci_block_size = fibBaseBits;

// The amount of blocks we need to encode the largest number if we use fibBase for storing all blocks
constexpr uint32 fibonacci_blocks_needed = (max_needed_fibonacci_index + fibonacci_block_size - 1) / fibonacci_block_size;

/* Build a precomputed table with enough size to hold all fibonacci numbers that are less
 * than the maximum number representable by our base type
 */

// Since we actually don't separate fibonacci(1) and fibonacci(2) (they are both equal to 1), we need to add 1,
// and because using n bits we can encode all numbers up to fibonacci(n+1), we add another 1 and get the next number after
// the largest number which we can encode using n bits (fibonacci_block_size bits in our case)
constexpr uint64 precalc_table_limit = FastFibonacci(fibonacci_block_size + 2);
constexpr std::array<fibBase, fibonacci_block_size> precalc_fibonacci_table = BuildPrecomputedTable<fibonacci_block_size>();

/* Lightweight class for storing a bitset of fixed size.
 *
 * This stores the bits in the order which is defined by the system (so for ex. on a little endian system this would store
 * little-endian numbers, specifically little-endian int16's by default, and with big-endian in each byte, since that is default).
 * Of course, this won't correspond to the encoding, but what matters is actually how we interpret the bytes, and because
 * we read and write them all in the same way, this should be perfectly fine. As an addition, writing and reading this way
 * lets us write cleaner code.
 *
 * It's easy to change this to a template class, but currently there's no need, cause we have a constant
 * maximum code length, and the overhead of a more complex class is bigger than just storing the extra bytes
 */
struct LightBitset {
	protected:
		// the max length is the index of the maximum fibonacci number we can use, rounded up to the block size
		static constexpr uint32 BLOCKS = fibonacci_blocks_needed;
		// the size of the block, it makes sense to use fibBase size for it in order for it to properly fit into the whole project
		// also this means that we can copy bits by blocks without even having to worry about overlapping or whatever (when constructing the fib encoding)
		static constexpr uint32 BLOCKSIZE = fibonacci_block_size;
		// For less time consuming operations (because otherwise there's a lot of useless division if we don't Ofast)
		static constexpr uint32 BLOCKSIZE_SHIFT = gcem::round(gcem::log2(BLOCKSIZE));
		static constexpr uint32 BLOCKSIZE_MASK = BLOCKSIZE - 1;
	public:
		fibBase bits_[BLOCKS];
		uint8 max_index_;

		constexpr LightBitset() {
			max_index_ = 0;
			Clear();
		}

		~LightBitset() = default;

		/* Basic functions for manipulating instances of this class */

		constexpr void Set(uint8 index) {
			bits_[(index >> BLOCKSIZE_SHIFT)] |= (1 << (index & BLOCKSIZE_MASK));
			if (index + 1 > max_index_)
				max_index_ = index + 1;
		}

		constexpr void SetBlock(uint8 block_index, fibBase block) {
			bits_[block_index] = block;
		}

		constexpr bool Get(uint8 index) const {
			return (bits_[(index >> BLOCKSIZE_SHIFT)] &  (1 << (index & BLOCKSIZE_MASK))) != 0;
		}

		constexpr fibBase GetBlock(uint8 block_index) {
			return bits_[block_index];
		}

		constexpr void Clear() {
			for (uint32 i = 0; i < BLOCKS; i++)
				bits_[i] = 0;
		}

		/* Function for merging this bitset with another one
		 *
		 * Pretty much just sets all the bits present in the other bitset
		 * Merges by block, not by bit. This is correct even if the last block in the other bitset
		 * isn't completely filled, because in that case the values will be 0 anyways
		 *
		 * destIndex specifies where to start putting the blocks from other into the current bitset
		 */
		constexpr void MergeWith(const LightBitset &other, const uint8 destIndex) {
			for (uint32 i = 0; i < ((other.max_index_ + BLOCKSIZE_MASK - 1) >> BLOCKSIZE_SHIFT); i++) {
				#ifdef DEBUG
				if (i + destIndex > BLOCKS) {
					puts("INVALID MERGING OF BLOCKS");
					exit(-1);
				}
				#endif
				bits_[i + destIndex] = (bits_[i + destIndex] | other.bits_[i]);
			}
			if (other.max_index_ + destIndex * BLOCKSIZE > max_index_)
				max_index_ = other.max_index_ + destIndex * BLOCKSIZE;
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
	int32 last = fibonacci_block_size;
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

// We have our one main table which stores all codes for numbers which can be written using fibBaseBits in fibonacci encoding,
// so we need to store a conversion table for the rest of the blocks, which we will use to quickly calculate correct shifts
// and encodings for the numbers in that range.
constexpr uint32 extra_fibonacci_blocks_needed = fibonacci_blocks_needed - 1;

/* Function which generates a "shift table" which we used for fast fibonacci encoding
 * The size of the table is m*n, where m is the number of extra blocks needed (number of total blocks - 1)
 * and m is the maximum number we can encode using fibonacci_block_size bits with fibonacci encoding.
 * The only data we actually store here is the range in which any number shifted by k*block_size "fibonacci bits"
 * would result in the number i (so pretty much [k][i] is a pair).
 */
template<size_t n, size_t m>
constexpr std::array<
	std::array<std::pair<uint64, uint64>, n>,
	m> MakeFibonacciShiftTable() {
		std::array<
			std::array<std::pair<uint64, uint64>, n>,
			m> result = {};
		// Create a table of mins and maxs for every shift (k*fibonacci_block_size)
		for (uint32 k = 0; k < m; k++) {
			// The fibonacci number which was the largest in the previous block
			uint64 previous_maximum = FastFibonacci((k + 1) * fibonacci_block_size + 1);
			// The fibonacci number which is the bare minimum for the current block
			uint64 current_minimum = FastFibonacci((k + 1) * fibonacci_block_size + 2);

			// We can't get a 0 because before shifting that the number actually needs to be shifted, i.e.
			// it is bigger than current_minimum. From here on the maximum value for current i is calculated
			// as the minimum + current_minimum - 1 if the result of the shift doesn't have the first bit set
			// (which would mean that we can take current_minimum - 1 extra numbers, starting with
			// fibonacci((k + 1) * fib_block_size - 1)), so the previous largest number and then all the ones lower than it
			// (which gives us the maximum value current_minimum - 1), however if it DOES have the first bit set,
			// then we add previous_maximum + 1, which corresponds to us not being to make current_minimum - 1 from the previous
			// bits, since we can't take the number previous_maximum, so we can actually make only previous_maximum - 1
			result[k][1].first = current_minimum;
			result[k][1].second = current_minimum + previous_maximum - 1;

			for (uint32 i = 2; i < n; i++) {
				result[k][i].first = result[k][i - 1].second + 1;
				result[k][i].second = result[k][i].first + (precalc_code_table[i].Get(0) ? previous_maximum : current_minimum) - 1;
			}
		}
		return result;
	}

// Table which stores all the necessary information for calculating the encoding of a number by shifting it's fibonacci
// representation. As it turns out, there's only use in storing the minimum and maximum number for which the result
// of the k*fibBaseBits'th shift is equal to the resulting number. As, during the calculation of the shift, we can get
// an inconsistency in the interval [0, 1], we first round and check if the number if it fits, and if it doesn't - subtract 1
// from it and get the actual proper representation.
constexpr std::array<
    std::array<std::pair<uint64, uint64>, precalc_table_limit>,
    extra_fibonacci_blocks_needed> fibonacci_shift_table =
    	MakeFibonacciShiftTable<precalc_table_limit, extra_fibonacci_blocks_needed>();

template <size_t n>
constexpr std::array<long double, n> GenerateRequiredPhiShifts() {
	std::array<long double, n> result = {};
	for (uint32 i = 0; i < n; i++) {
		result[i] = gcem::pow(phi, (i + 1) * fibonacci_block_size);
	}
	return result;
}

constexpr std::array<long double, extra_fibonacci_blocks_needed> golden_ratio_shift_powers =
	GenerateRequiredPhiShifts<extra_fibonacci_blocks_needed>();

constexpr LightBitset EncodeNumber(uint64 n) {
	if (n < precalc_table_limit)
		return precalc_code_table[n];
	LightBitset result;
	for (int32 i = extra_fibonacci_blocks_needed - 1; i >= 0; i--) {
		if (n >= fibonacci_shift_table[i][1].first) {
			uint64 current = gcem::round(static_cast<long double>(n) / golden_ratio_shift_powers[i]);
			if (n < fibonacci_shift_table[i][current].first) {
				current--;
			}
			#ifdef DEBUG
			// simple check for debuggind, second if case is for the numbers which are so close to 2^64 that the right side has overflow,
			// which wouldn't matter in the case of actual encoding, but could mess up the if logic here
			if (!(fibonacci_shift_table[i][current].first <= n && fibonacci_shift_table[i][current].second >= n) &&
					!(fibonacci_shift_table[i][current].second < fibonacci_shift_table[i][current].first)) {
				printf("INVALID ENCODING OF NUMBER %llu\n", n);
				exit(-1);
			}
			#endif
			n -= fibonacci_shift_table[i][current].first;
			result.MergeWith(precalc_code_table[current], i + 1);
		}
	}

	#ifdef DEBUG
	if (n >= precalc_table_limit) {
		printf("INVALID NUMBER %llu AS A RESULT OF DOING SHIFTS, CAN'T ENCODE\n", n);
		exit(-1);
	}
	#endif

	if (n > 0) {
		result.MergeWith(precalc_code_table[n], 0);
	}

	return result;
}

int main() {
	LightBitset encoded = EncodeNumber(2178309);
	for (uint32 j = 0; j < encoded.max_index_; j++) {
		std::cout << (encoded.Get(j) ? '1' : '0');
	}
	std::cout << '1' << std::endl;

	std::cout << fibonacci_block_size << " " << fibonacci_blocks_needed << std::endl;
	std::cout << precalc_table_limit << std::endl;

	for (uint32 i = 0; i < fibonacci_block_size; i++) {
		std::cout << i << ": " << precalc_fibonacci_table[i] << std::endl;
	}
//	while (true) {
//		int i;
//		std::cin >> i;
//		std::cout << MaxFittingFibonacci(i) << " " << FibonacciEncodingLength(i) << std::endl;
//	}

	for (uint32 i = 0; i < 150; i++) {
		LightBitset encoded = precalc_code_table[i];
		std::cout << i << " ";
		for (uint32 j = 0; j < encoded.max_index_; j++) {
			std::cout << (encoded.Get(j) ? '1' : '0');
		}
		std::cout << '1' << std::endl;
	}

	return 0;
}
