// SPDX-FileCopyrightText: Copyright (c) 2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DMRG_PERMUTATION_PARITY_H
#define DMRG_PERMUTATION_PARITY_H

#include <PsimagLite/Vector.h>
#include <vector>

namespace Dmrg {

inline int parityOfPermutation(const std::vector<SizeType>& permutation) noexcept
{
	int parity = 1;
	for (SizeType i = 0; i < permutation.size(); ++i)
		for (SizeType j = i + 1; j < permutation.size(); ++j)
			if (permutation[i] > permutation[j])
				parity = -parity;

	return parity;
}

} // namespace Dmrg

#endif // DMRG_PERMUTATION_PARITY_H
