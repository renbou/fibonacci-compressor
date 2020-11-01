#include <iostream>
#include "fibonacci_encoding.hpp"
#include "pretty_types.h"

int main() {
//	const FibonacciBitset &encoded = EncodeNumber(25799813418236517);
//	std::cout << 25799813418236517 << " ";
//	for (uint32 j = 0; j < encoded.max_index_; j++) {
//		std::cout << (GetFibonacciBitsetBit(encoded, j) ? '1' : '0');
//	}
//	std::cout << '1' << std::endl;

	srandom(time(nullptr));
	clock_t start = clock();
	for (uint32 i = 0; i < 100; i++) {
		uint64 number = random() * random();
		const FibonacciBitset &encoded = EncodeNumber(number);
		std::cout << number << " ";
		for (uint32 j = 0; j < encoded.max_index_; j++) {
			std::cout << ((encoded.bits_[j >> fibonacci_blocksize_shift] & (1 << (j & fibonacci_blocksize_mask))) ? '1' : '0');
		}
		std::cout << '1' << std::endl;
		if ((encoded.bits_[0] & (1 << 12)) == 1 && (encoded.bits_[2] & (1 << 8)) == 0) {
			srandom(time(nullptr));
		}
	}
	clock_t stop = clock();
	double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
	printf("\nTime elapsed: %.5f\n", elapsed);

//
//	for (uint32 i = 0; i < 150; i++) {
//		const FibonacciBitset &encoded = EncodeNumber(i);
//		std::cout << i << " ";
//		for (uint32 j = 0; j < encoded.max_index_; j++) {
//			std::cout << (GetFibonacciBitsetBit(encoded, j) ? '1' : '0');
//		}
//		std::cout << '1' << std::endl;
//	}

	return 0;
}
