// Created by Artem Mikheev on 1/11/20.
// Apache License 2.0
#ifndef PHICC_FIBONACCI_ENCODING_H
#define PHICC_FIBONACCI_ENCODING_H

#include <array>
#include <gcem.hpp>
#include "pretty_types.h"

/* * * * * * * * * * * * * * *
 * Basic constants and types *
 * * * * * * * * * * * * * * */

/* fibBase is the type which we will be using to store 1 block of the fibonacci encoded number */
typedef uint16 fibBase;
/* how many bits we can actually store, needed to understand what's the max number we can encode with 1 block */
constexpr uint32 fibBaseBits = sizeof(fibBase) * 8;
/* The size of every block of our encoding, this is just for simpler code reading
 * and should always be equal to fibBaseBits */
constexpr uint32 fibonacci_block_size = fibBaseBits;

/* phi and its reciprocal, needed for formulas */
constexpr long double phi = (1.0L + gcem::sqrt(5.0L)) / 2.0L;
constexpr long double Phi = 1.0L - phi;

/* Struct for storing one encoded block, in reality this will only be used in the precomputed table */
struct FibonacciBlock {
	fibBase bits_;
	uint8 max_index_;

	constexpr FibonacciBlock()
	: bits_(0), max_index_(0) {}
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Various functions for helping with fibonacci coding *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* This file implements all of the functions it defines in order for them to be constexpr-valid,
 * so sadly this isn't insanely beautiful
 */

/* Calculates i'th fibonacci number very quickly using Binet's formula */
constexpr uint64 FastFibonacci(const uint64 i) {
	return (gcem::pow(phi, i) - gcem::pow(Phi, i)) / gcem::sqrt(5.0L);
}

/* Function that returns the index of the largest fibonacci number that's less than n
 * (also the length of n encoded with fibonacci)
 *
 * This implementation can be derived using Binet's representation and the fact that (1 - phi)^n -> 0 as n -> inf
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
constexpr uint32 MaxFittingFibonacciIndex(uint64 n) {
	if (n < 20)
		return small_number_fitter[n];
	constexpr long double sqrt5 = gcem::sqrt(5.0L);
	constexpr long double denominator = gcem::log(1.0L + sqrt5) - gcem::log(2.0L);
	// add some epsilon to escape invalid answers for fibonacci numbers
	// higher than the ones stored in the small fitter
	return gcem::log(sqrt5 * static_cast<long double>(n)) / denominator + 0.001L;
}

/* Calculates the length of the fibonacci encoding of decimal n
 *
 * Pretty much just a trunk for the last function, but handles 0 differently
 * because the encoding of 0 is nothing (actually just the '1' which is placed after every code)
 */
constexpr inline uint32 FibonacciEncodingLength(uint64 n) {
	if (n != 0) [[likely]] {
		return MaxFittingFibonacciIndex(n) - 1;
	}
	return 0;
}

/* Build precomputed at fibonacci table of size n
 * (so the first n fibonacci numbers not counting 1 cause its the same as 2)
 *
 * We later use this in order to build the base encoding table for small values
 * As well as for simply converting small numbers to fibonacci encoding
 *
 * If the n passed is larger than the maximum encodable by T, this will just be wrong
 */
template <typename T, size_t n>
constexpr std::array<T, n> BuildPrecomputedFibonacciTable() {
	std::array<T, n> result = {1};
	if (n > 1)
		result[1] = 2;
	for (T i = 2; i < n; i++)
		result[i] = result[i - 1] + result[i - 2];
	return result;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Stuff for encoding smaller numbers and creating the first precomputed tables for encoding stuff *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/* Length of the maximum code we can get, we only use up to uint64 so this will be enough for everything */
constexpr uint32 max_needed_fibonacci_index = MaxFittingFibonacciIndex(std::numeric_limits<uint64>::max());

/* The amount of blocks we need to encode the largest number using our block size */
constexpr uint32 max_needed_fibonacci_blocks = (max_needed_fibonacci_index + fibonacci_block_size - 1) / fibonacci_block_size;

/* The size of the precalculated table which we will use for encoding small numbers
 *
 * This is equivalent to the next number after the largest one we can encode using one block.
 * We +2 here due to the fact that we skipp fib(1) since it is equal to fib(2), and also because
 * we start our table from index 0, not 1.
 */
static constexpr uint32 max_encodable_by_one_block = FastFibonacci(fibonacci_block_size + 2);

/* Build the precalculated table with encodings for all the numbers which we can encode with one fibBase block */
constexpr std::array<fibBase, fibonacci_block_size> precalc_fibonacci_table =
	BuildPrecomputedFibonacciTable<fibBase, fibonacci_block_size>();

// For less time consuming operations (because otherwise there's a lot of useless division if we don't Ofast)
static constexpr uint32 fibonacci_blocksize_shift = gcem::round(gcem::log2(fibonacci_block_size));
static constexpr uint32 fibonacci_blocksize_mask = fibonacci_block_size - 1;

//\\-----------------------------------\\//

/* Encodes a small number into its fibonacci encoding
 *
 * N must be less than the maximum number encodable with fibBaseBits bits of fibonacci encoding,
 * otherwise we will get the wrong encoding
 */
constexpr FibonacciBlock EncodeSmallNumber(fibBase n) {
	FibonacciBlock result;
	int32 last = fibonacci_block_size;
	while (n > 0) {
		int32 left = 0, right = last;
		while (left < right) {
			int64 mid = (left + right) >> 1;
			if (precalc_fibonacci_table[mid] > n)
				right = mid;
			else
				left = mid + 1;
		}
		last = right;
		n -= precalc_fibonacci_table[right - 1];

		/* Set the required bit in the result and update max index */
		int32 index = right - 1;
		result.bits_ |= (1 << (index & fibonacci_blocksize_mask));
		if (index + 1 > result.max_index_)
			result.max_index_ = index + 1;
	}
	return result;
}

/* This generates a code table for small numbers which we can encode using just 1 fibBase block
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
constexpr std::array<FibonacciBlock, n> GenerateCodeTable() {
	std::array<FibonacciBlock, n> result = {};
	for (fibBase i = 0; i < n; i++) {
		result[i] = EncodeSmallNumber(i);
	}
	return result;
}

/* Build a table of all the codes we can encode with just 1 fibBase block */
static constexpr std::array<FibonacciBlock, max_encodable_by_one_block> precomputed_code_table =
	GenerateCodeTable<max_encodable_by_one_block>();


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Stuff for encoding numbers of any valid size and building all the required precomputed tables *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/* Struct for storing multiple encoded blocks, it would be quite slow to constantly copy this,
 * so actually we only store one global instance of this and always return the pointer to it.
 * IF YOU WANT TO ENCODE MULTIPLE DATA STREAMS AT ONE TIME THEN YOU MUST COPY THIS YOURSELF,
 * DO NOT RELY ON THIS LIBRARY RETURNING IT IN A SAFE MANNER TO USE INSIDE AN ARRAY OR WHATEVER
 */
struct FibonacciBitset {
	std::array<fibBase, max_needed_fibonacci_blocks> bits_;
	uint8 max_index_;

	constexpr FibonacciBitset()
	: bits_(), max_index_(0) {}
};

static FibonacciBitset global_fibonacci_bitset;

/* Helper function mostly for testing which gets the i'th bit in the bitset */
inline bool GetFibonacciBitsetBit(const FibonacciBitset &bitset, uint32 index) {
#ifdef DEBUG
	if (!((index >> fibonacci_blocksize_shift) < max_needed_fibonacci_blocks)) {
		printf("Accessing invalid index %u of bitset\n", index);
		exit(-1);
	}
#endif
	return bitset.bits_[index >> fibonacci_blocksize_shift] & (1 << (index & fibonacci_blocksize_mask));
}

/* Function for merging a bitset with the specified block, which is to be placed at the specified index */
inline void MergeFibonacciBitsetWithBlock(FibonacciBitset &bitset, const FibonacciBlock &block, uint32 index) {
	bitset.bits_[index] = block.bits_;
	if (index * fibonacci_block_size + block.max_index_ > bitset.max_index_) {
		bitset.max_index_ = index * fibonacci_block_size + block.max_index_;
	}
}

/* Function which generates a "shift table" which we use for fast fibonacci encoding
 *
 * The size of the table is shifts_n*n, where shifts_ is the number of extra blocks needed (number of total blocks - 1)
 * and n is the maximum number we can encode using fibonacci_block_size bits with fibonacci encoding.
 * The only data we actually store here is the range in which any number shifted by k*block_size "fibonacci bits"
 * would result in the number i (so pretty much [k][i] is a pair), which is needed because when we use our
 * so-called "fast fibonacci shift", we can lose due to the rounding up, and need to check if we actually got
 * the correct range. Of course we could calculate it every time using fibonacci(k)*n where k is the shift
 * (so for us its 16, 32, 48, 64, ...), but it would be quite slow to use quite slow maths every time for this.
 */
template<size_t n, size_t shifts_n>
constexpr std::array<
	std::array<std::pair<uint64, uint64>, n>,
	shifts_n> MakeFibonacciShiftTable()  {
	std::array<
		std::array<std::pair<uint64, uint64>, n>,
		shifts_n> result = {};
	// Create a table of mins and maxs for every shift (k*fibonacci_block_size)
	for (uint32 k = 0; k < shifts_n; k++) {
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
			result[k][i].second = result[k][i].first + ((precomputed_code_table[i].bits_ & 1) ? previous_maximum : current_minimum) - 1;
		}
	}
	return result;
}

/* Create a simple table of the required fibonacci-shifts
 *
 * Which are all actually phi^some_power where some_power are actually the multiples of our block size,
 * so 16, 48, 64, ... We don't save the 0'th shift because it's quite stupid, since phi^0 = 1, and also
 * we never use it, because we have a prebuilt table for all the number which are encoded using 1 block.
 */
template <size_t shifts_n>
constexpr std::array<long double, shifts_n> GenerateRequiredPhiShifts()  {
	std::array<long double, shifts_n> result = {};
	for (uint32 i = 0; i < shifts_n; i++) {
		result[i] = gcem::pow(phi, (i + 1) * fibonacci_block_size);
	}
	return result;
}

//\\-----------------------------------\\//

/* How many extra shifts we would need in order to encode all the numbers
 *
 * This is fibonacci_blocks_needed - 1 because we have our precomputed table for the first block
 */
static constexpr uint32 extra_fibonacci_blocks_needed = max_needed_fibonacci_blocks - 1;

/* Table which stores all the necessary information for calculating the encoding of a number by shifting it's fibonacci
 * representation. As it turns out, there's only use in storing the minimum and maximum number for which the result
 * of the k*fibBaseBits'th shift is equal to the resulting number. As, during the calculation of the shift, we can get
 * an inconsistency in the interval [0, 1], we first round and check if the number if it fits, and if it doesn't - subtract 1
 * from it and get the actual proper representation.
 */
static constexpr std::array<
	std::array<std::pair<uint64, uint64>, max_encodable_by_one_block>,
	extra_fibonacci_blocks_needed> fibonacci_shift_table =
	MakeFibonacciShiftTable<max_encodable_by_one_block, extra_fibonacci_blocks_needed>();

/* Table holding all of the required powers of phi for fast shifting
 */
static constexpr std::array<long double, extra_fibonacci_blocks_needed> phi_shift_powers =
	GenerateRequiredPhiShifts<extra_fibonacci_blocks_needed>();

//\\-----------------------------------\\//

/* Returns a pointer to the encoding of the number n
 *
 * Current implementations agrees upon the fact that the returned pointer shall always be the same and const,
 * because we use one single instance of this struct for encoding all numbers, in order to save ourselves from
 * massive amounts of useless C++ copy constructors, just to make the compressor faster.
 */
const FibonacciBitset& EncodeNumber(uint64 n)  {
#ifdef DEBUG
	uint64 original = n;
#endif
	// Reset the bitset pretty much
	global_fibonacci_bitset.max_index_ = 0;
	if (n < max_encodable_by_one_block) {
		global_fibonacci_bitset.bits_[0] = precomputed_code_table[n].bits_;
		global_fibonacci_bitset.max_index_ = precomputed_code_table[n].max_index_;
		return global_fibonacci_bitset;
	}

	for (int32 i = extra_fibonacci_blocks_needed - 1; i >= 0; i--) {
		if (n >= fibonacci_shift_table[i][1].first) {
			uint64 current = gcem::round(static_cast<long double>(n) / phi_shift_powers[i]);
			// We might get a round up into the wrong direction, resulting in either the next number OR the limit,
			// which we must check first in order to avoid out-of-bounds.
			if (current == max_encodable_by_one_block || n < fibonacci_shift_table[i][current].first ) {
				current--;
			}
#ifdef DEBUG
			// simple check for debugging, second if case is for the numbers which are so close to 2^64 that the right side has overflow,
			// which wouldn't matter in the case of actual encoding, but could mess up the if logic here
			if (!((fibonacci_shift_table[i][current].first <= n && fibonacci_shift_table[i][current].second >= n) ||
						   (fibonacci_shift_table[i][current].second < fibonacci_shift_table[i][current].first))) {
				printf("INVALID NUMBER %llu AS A RESULT OF DOING SHIFTS TO, CAN'T ENCODE %llu\n", n, original);
				exit(-1);
			}
#endif
			n -= fibonacci_shift_table[i][current].first;
			MergeFibonacciBitsetWithBlock(global_fibonacci_bitset, precomputed_code_table[current], i + 1);
		}
	}

#ifdef DEBUG
	if (!(n < max_encodable_by_one_block)) {
		printf("INVALID NUMBER %llu AS A RESULT OF DOING SHIFTS TO, CAN'T ENCODE %llu\n", n, original);
		exit(-1);
	}
#endif

	if (n > 0) {
		global_fibonacci_bitset.bits_[0] = precomputed_code_table[n].bits_;
		MergeFibonacciBitsetWithBlock(global_fibonacci_bitset, precomputed_code_table[n], 0);
	}

	return global_fibonacci_bitset;
}

#endif //PHICC_FIBONACCI_ENCODING_H
