// SPDX-FileCopyrightText: Copyright (c) 2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PermutationParity.h"
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <numeric>
#include <random>
#include <vector>

TEST_CASE("Permutation parity handles identity and basic cycles", "[PermutationParity]")
{
	CHECK(Dmrg::parityOfPermutation(std::vector<SizeType> {}) == 1);
	CHECK(Dmrg::parityOfPermutation(std::vector<SizeType> { 0 }) == 1);
	CHECK(Dmrg::parityOfPermutation(std::vector<SizeType> { 0, 1, 2, 3 }) == 1);
	CHECK(Dmrg::parityOfPermutation(std::vector<SizeType> { 1, 0 }) == -1);
	CHECK(Dmrg::parityOfPermutation(std::vector<SizeType> { 1, 2, 0 }) == 1);
	CHECK(Dmrg::parityOfPermutation(std::vector<SizeType> { 1, 0, 3, 2 }) == 1);
}

TEST_CASE("Permutation parity handles reverse-order permutations", "[PermutationParity]")
{
	for (SizeType size = 0; size <= 16; ++size) {
		std::vector<SizeType> permutation(size);
		std::iota(permutation.begin(), permutation.end(), 0);
		std::reverse(permutation.begin(), permutation.end());

		const SizeType inversions = (size < 2) ? 0 : size * (size - 1) / 2;
		const int      expected   = (inversions % 2 == 0) ? 1 : -1;
		CHECK(Dmrg::parityOfPermutation(permutation) == expected);
	}
}

TEST_CASE("Permutation parity is not limited by an integer bit mask", "[PermutationParity]")
{
	std::vector<SizeType> permutation(80);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::swap(permutation[64], permutation[65]);

	CHECK(Dmrg::parityOfPermutation(permutation) == -1);
}

TEST_CASE("Permutation parity follows randomized transposition sequences", "[PermutationParity]")
{
	std::mt19937                  rng(1729);
	const std::array<SizeType, 7> sizes { 2, 3, 4, 8, 17, 65, 96 };

	for (const SizeType size : sizes) {
		std::uniform_int_distribution<SizeType> index(0, size - 1);
		for (SizeType trial = 0; trial < 25; ++trial) {
			std::vector<SizeType> permutation(size);
			std::iota(permutation.begin(), permutation.end(), 0);
			int expected = 1;

			const SizeType swaps = 73 + trial;
			for (SizeType step = 0; step < swaps; ++step) {
				const SizeType first  = index(rng);
				SizeType       second = index(rng);
				while (second == first)
					second = index(rng);
				std::swap(permutation[first], permutation[second]);
				expected = -expected;
			}

			CHECK(Dmrg::parityOfPermutation(permutation) == expected);

			std::vector<SizeType> inverse(size);
			for (SizeType i = 0; i < size; ++i)
				inverse[permutation[i]] = i;
			CHECK(Dmrg::parityOfPermutation(inverse) == expected);
		}
	}
}
