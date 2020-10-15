#include <iostream>
#include <vector>
#include <intx/intx.hpp>
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

typedef uint64 fibBase;
typedef uint32 fibBaseHalf;
#define FIB_BASE_BITS 64
static_assert(FIB_BASE_BITS == sizeof(fibBase) * 8, "preprocessor-time FIB_BASE_BITS isn't equal to compile-time fibBaseBits");
static_assert(sizeof(fibBaseHalf) * 2 == sizeof(fibBase), "fibBaseHalf must be a type with 2 times less bits than fibBase");
constexpr fibBase MAX_FIBBASE = std::numeric_limits<fibBase>::max();
constexpr fibBase MAX_FIBBASEHALF = std::numeric_limits<fibBaseHalf>::max();
// I couldn't figure out how to get the preprocessor to use compile-time constexpr values (makes sense I guess)
// So in order for the compilation to work properly, if you change this to intx::uint128 or more, you should
// Set FIB_BASE_BITS accordingly so that the appropriate constexpr functions will be used


// constexpr function for counting the logarithm at compile time, taken from some old code of mine
constexpr uint64 constexprLog2(uint64 n) {
	return n <= 1 ? 0 : 1 + constexprLog2(n / 2);
}

// number of bits in the base, we add one to count for the 0'th
constexpr uint64 fibBaseBits = sizeof(fibBase) * 8;
constexpr uint64 fibBaseBitsPower = constexprLog2(fibBaseBits) + 1;

/* This matrix type will always be the same size, only the data type might change
 *
 * Actually only used for fast fibonacci number calculation in log(n) using binary matrix exponentiation
 */
typedef std::array<std::array<fibBase, 2>, 2> matrix;

/* Hacky constexpr function that we will use everywhere for multiplication
 * Needed because operator * for classes isn't constexpr, but of course the multiplication itself is
 */
constexpr fibBase constexpr_multiply(const fibBase &a, const fibBase &b) {
#if FIB_BASE_BITS <= 64
	return a * b;
#elif FIB_BASE_BITS == 128
	return intx::constexpr_mul(a, b);
#else
	return intx::constexpr_mul<FIB_BASE_BITS>(a, b);
#endif
}

/* Hacky constexpr functions that will get the lower and upper halves of the number, no matter the size
 * Needed because operator << for classes isn't constexpr
 */
constexpr fibBaseHalf constexpr_lo(const fibBase &a) {
#if FIB_BASE_BITS <= 64
	return a & MAX_FIBBASEHALF;
#else
	return a.lo;
#endif
}

constexpr fibBaseHalf constexpr_hi(const fibBase &a) {
#if FIB_BASE_BITS <= 64
	return a >> (fibBaseBits / 2);
#else
	return a.hi;
#endif
}

// Simple function which multiplies two 2-d matrices, trying to evaluate this value at compile time
constexpr matrix MatrixMultiply(const matrix &a, const matrix &b) {
	return matrix {
		{{constexpr_multiply(a[0][0], b[0][0]) + constexpr_multiply(a[0][1], b[1][0]),
			constexpr_multiply(a[0][0], b[0][1]) + constexpr_multiply(a[0][1], b[1][1])},
   		{constexpr_multiply(a[1][0], b[0][0]) + constexpr_multiply(a[1][1], b[1][0]),
	  		constexpr_multiply(a[1][0], b[0][1]) + constexpr_multiply(a[1][1], b[1][1])}}
	};};

//#if FIB_BASE_BITS <= 64
//constexpr matrix MatrixMultiply(const matrix &a, const matrix &b) {
//	return matrix {{{a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]},
//					   {a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]}}};
//}
//#elif FIB_BASE_BITS == 128
////constexpr matrix MatrixMultiply(const matrix &a, const matrix &b) {
////	return matrix {
////		{{intx::constexpr_mul(a[0][0], b[0][0]) + intx::constexpr_mul(a[0][1], b[1][0]),
////			intx::constexpr_mul(a[0][0], b[0][1]) + intx::constexpr_mul(a[0][1], b[1][1])},
////   		 {intx::constexpr_mul(a[1][0], b[0][0]) + intx::constexpr_mul(a[1][1], b[1][0]),
////		 	intx::constexpr_mul(a[1][0], b[0][1]) + intx::constexpr_mul(a[1][1], b[1][1])}}};
////}
//#else
//constexpr matrix MatrixMultiply(const matrix &a, const matrix &b) {
//	return matrix {
//		{{intx::constexpr_mul<FIB_BASE_BITS>(a[0][0], b[0][0]) + intx::constexpr_mul<FIB_BASE_BITS>(a[0][1], b[1][0]),
//			intx::constexpr_mul<FIB_BASE_BITS>(a[0][0], b[0][1]) + intx::constexpr_mul<FIB_BASE_BITS>(a[0][1], b[1][1])},
//   		 {intx::constexpr_mul<FIB_BASE_BITS>(a[1][0], b[0][0]) + intx::constexpr_mul<FIB_BASE_BITS>(a[1][1], b[1][0]),
//		 	intx::constexpr_mul<FIB_BASE_BITS>(a[1][0], b[0][1]) + intx::constexpr_mul<FIB_BASE_BITS>(a[1][1], b[1][1])}}};
//}
//#endif

// Function which calculates the first n matrix powers in the form 2^i
// (so base_matrix^1, base_matrix^2, base_matrix^4, ...)
template <size_t n>
constexpr std::array<matrix, n> CalculateMatrixPowers(const matrix& base_matrix) {
	std::array<matrix, n> result = {base_matrix};
	// Use matrix multiplication to calculate the matrices, storing them by rvalue reference
	for (uint32 i = 1; i < n; i++) {
		result[i] = MatrixMultiply(result[i-1], result[i-1]);
	}
	return result;
}

// The base fibonacci matrix, this is used for all later calculations
static constexpr matrix fibonacci_matrix = {{{1, 1}, {1, 0}}};
constexpr std::array<matrix, fibBaseBitsPower> first_fibmatrix_powers = CalculateMatrixPowers<fibBaseBitsPower>(fibonacci_matrix);
// Simply get the actual values from the matrices, we save the matrices
// In case we later need them for fast fibonacci calculation
constexpr std::array<fibBase, fibBaseBitsPower> first_fibonacci_powers = [](){
	std::array<fibBase, fibBaseBitsPower> result = {first_fibmatrix_powers[0][0][1]};
	for (uint32 i = 0; i < fibBaseBitsPower; i++) {
		result[i] = first_fibmatrix_powers[i][0][1];
	}
	return result;
}();

constexpr bool WillOverflow(fibBase a, fibBase b) {
	if (a > b) {
		fibBase tmp = std::move(a);
		a = std::move(b);
		b = std::move(tmp);
	}
	if (a > MAX_FIBBASEHALF)
		return true;
	fibBaseHalf c = constexpr_hi(b);
	if (constexpr_multiply(a, c) > MAX_FIBBASEHALF)
		return true;
	return false;
}

constexpr bool MatrixRowMultWillOverflow(fibBase row1, fibBase row2, fibBase col1, fibBase col2) {
	return WillOverflow(row1, col1) || WillOverflow(row2, col2) ||
		(constexpr_multiply(row1, col1) + constexpr_multiply(row2, col2)) < constexpr_multiply(row1, col1);
}

/* Function that returns the index of the maximum fibonacci number which is less than n
 * We use this in order to construct a basic table of fibonacci numbers
 *
 * In order to search for the fibonacci number, we greedily take the largest one we can
 * on every step, what this means is we try to multiply the previously gotten result (on the last step)
 * by a possible matrix, getting the prev_power+i'th fibonacci number. If the fibonacci number overflows
 * or is greater than what we need, we continue searching.
 */
constexpr uint64 MaxFittingFibonacci(fibBase n) {
	// initialize result as identity matrix for the first round
	matrix result = {{{1, 0}, {0, 1}}};
	int64 i = fibBaseBitsPower, result_i = 0;
	// variable to mark that the next fibonacci number has been overflown
	bool end_is_near = false;
	while (i > 0 && !end_is_near) {
		// search for the next item using binary search, this is just in order to optimize out lots of the matrix multiplicationss
		int64 left = 0, right = i - 1;
		while (left <= right) {
			int64 mid = (left + right) >> 1;

			// If we have overflown, then there's no point in doing anything else
			if (MatrixRowMultWillOverflow(result[0][0], result[0][1], first_fibmatrix_powers[mid][0][1], first_fibmatrix_powers[mid][1][1]))
				right = mid - 1;
			// separate this case for optimization, we don't need to multiply if we overflow anyway
			else if (MatrixMultiply(result, first_fibmatrix_powers[mid])[0][1] > n)
				right = mid - 1;
			else
				left = mid + 1;
		}

		// found a suitable item
		i = right;
		if (right >= 0) {
			// we have to check that if the next fibonacci number (in matrix it's at index 0,0) will be overflown during calculations,
			// then this is the end, and after this we have to stop
			end_is_near = MatrixRowMultWillOverflow(result[0][0], result[0][1], first_fibmatrix_powers[i][0][0], first_fibmatrix_powers[i][1][0]);

			result = MatrixMultiply(result, first_fibmatrix_powers[i]);
			result_i += (1 << i);
		}
	}
	return result_i;
}

uint64 precount_table_size = MaxFittingFibonacci(MAX_FIBBASE);

int main() {
//	for (uint32 i = 0; i < fibBaseBitsPower; i++) {
//		std::cout << intx::to_string(first_fibmatrix_powers[i][0][0]) << " " << intx::to_string(first_fibmatrix_powers[i][0][1]) << "\n" <<
//			intx::to_string(first_fibmatrix_powers[i][1][0]) << " " << intx::to_string(first_fibmatrix_powers[i][1][1]) << std::endl;
//	}
	std::cout << MaxFittingFibonacci(MAX_FIBBASE) << std::endl;
	for (uint32 i = 0; i < fibBaseBitsPower; i++) {
		std::cout << first_fibmatrix_powers[i][0][0] << " " << first_fibmatrix_powers[i][0][1] << "\n" <<
			first_fibmatrix_powers[i][1][0] << " " << first_fibmatrix_powers[i][1][1] << std::endl;
	}
	return 0;
}
