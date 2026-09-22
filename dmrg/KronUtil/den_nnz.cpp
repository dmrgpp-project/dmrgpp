#include "util.h"

template <typename ComplexOrRealType> int denNnz(const PsimagLite::Matrix<ComplexOrRealType>& a_)
{
	const int nrow_a = a_.n_row();
	const int ncol_a = a_.n_col();
	/*
	 * -------------------------
	 * return number of nonzeros
	 * matrix A in dense storage format
	 * -------------------------
	 */
	const bool use_estimate = true;
	if (use_estimate) {
		return (nrow_a * ncol_a);
	};

	const ComplexOrRealType dzero = 0;
	int                     nnz_a = 0;

	int ja = 0;

	/*
	 * -----------------------------
	 * Note that  nnz_A is a reduction variable
	 * -----------------------------
	 */
	for (ja = 0; ja < ncol_a; ja++) {
		int ia = 0;
		for (ia = 0; ia < nrow_a; ia++) {
			int is_zero = (a_(ia, ja) == dzero);
			nnz_a       = (is_zero) ? nnz_a : (nnz_a + 1);
		};
	};

	return (nnz_a);
}
#undef A
