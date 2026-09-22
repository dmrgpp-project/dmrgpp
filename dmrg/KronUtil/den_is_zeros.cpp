#include "util.h"

template <typename ComplexOrRealType>
bool denIsZeros(const PsimagLite::Matrix<ComplexOrRealType>& a_)
{
	const int nrow_a = a_.n_row();
	const int ncol_a = a_.n_col();
	/*
	 * -------------------------
	 * return whether A is the zero matrix
	 * matrix A in dense storage format
	 * -------------------------
	 */

	const ComplexOrRealType zero = 0;

	int ja = 0;
	for (ja = 0; ja < ncol_a; ja++) {
		int ia = 0;
		for (ia = 0; ia < nrow_a; ia++) {
			ComplexOrRealType aij = a_(ia, ja);
			if (aij != zero) {
				return (false);
			};
		};
	};

	/*
	 * -----------------
	 * passed all checks
	 * -----------------
	 */
	return (true);
}
#undef A
