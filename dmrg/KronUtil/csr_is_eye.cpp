#include "util.h"

template <typename ComplexOrRealType>
bool csrIsEye(const PsimagLite::CrsMatrix<ComplexOrRealType>& a)
{
	// ----------------------------------------------------
	// check whether a sparse matrix is the identity matrix
	// ----------------------------------------------------

	const int nrow_a = a.rows();
	const int ncol_a = a.cols();

	bool is_eye = (nrow_a == ncol_a);

	if (!is_eye) {
		return (false);
	};

	if ((nrow_a <= 0) || (ncol_a <= 0)) {
		return (false);
	};

	for (int ia = 0; ia < nrow_a; ia++) {
		int istart = a.getRowPtr(ia);
		int iend   = a.getRowPtr(ia + 1);

		bool has_diagonal = false;
		;
		for (int k = istart; k < iend; k++) {
			int               ja  = a.getCol(k);
			ComplexOrRealType aij = a.getValue(k);
			ComplexOrRealType eij = (ia == ja) ? 1 : 0;
			if (ia == ja) {
				has_diagonal = true;
			};

			is_eye = (aij == eij);
			if (!is_eye) {
				return (false);
			};
		};
		if (!has_diagonal) {
			return (false);
		};
	};

	return (true);
}
