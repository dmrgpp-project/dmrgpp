#include "util.h"

template <typename ComplexOrRealType>
void denSubmatrix(const int                                    nrow_A,
                  const int                                    ncol_A,
                  const PsimagLite::Matrix<ComplexOrRealType>& a_,
                  const int                                    nrindex,
                  const int                                    ncindex,
                  const PsimagLite::Vector<int>::Type&         rindex,
                  const PsimagLite::Vector<int>::Type&         cindex,
                  PsimagLite::Matrix<ComplexOrRealType>&       c_)
{
	/*
	 * -------------------------------------
	 * extract submatrix from a dense matrix
	 * -------------------------------------
	 */
	const int nrow_c = nrindex;
	const int ncol_c = ncindex;

	int ic = 0;
	int jc = 0;

	/*
	 * ------------------------
	 * check rindex[], cindex[]
	 * ------------------------
	 */
#ifndef NDEBUG
	for (ic = 0; ic < nrindex; ic++) {
		int ia   = rindex[ic];
		int isok = (0 <= ia) && (ia < nrow_A);
		assert(isok);
	};

	for (jc = 0; jc < ncindex; jc++) {
		int ja   = cindex[jc];
		int isok = (0 <= ja) && (ja < ncol_A);
		assert(isok);
	};
#endif

	for (jc = 0; jc < ncol_c; jc++) {
		for (ic = 0; ic < nrow_c; ic++) {
			int ia = rindex[ic];
			int ja = cindex[jc];

			c_(ic, jc) = a_(ia, ja);
		};
	};
}
