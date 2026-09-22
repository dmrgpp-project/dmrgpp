#include "util.h"

template <typename ComplexOrRealType>
void estimateKronCost(const int                                                nrow_A,
                      const int                                                ncol_A,
                      const int                                                nnz_A_in,
                      const int                                                nrow_B,
                      const int                                                ncol_B,
                      const int                                                nnz_B_in,
                      ComplexOrRealType*                                       p_kron_nnz,
                      ComplexOrRealType*                                       p_kron_flops,
                      int*                                                     p_imethod,
                      const typename PsimagLite::Real<ComplexOrRealType>::Type denseFlopDiscount)
{
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;

	int imethod = 1;

	assert(p_kron_nnz != 0);
	assert(p_kron_flops != 0);
	assert(p_imethod != 0);

	RealType nnz_a      = static_cast<RealType>(nnz_A_in);
	RealType nnz_b      = static_cast<RealType>(nnz_B_in);
	bool     is_dense_a = (nnz_A_in == nrow_A * ncol_A);
	bool     is_dense_b = (nnz_B_in == nrow_B * ncol_B);

	/*
	 * ------------------------------------
	 * assume dense matrix operations are
	 * faster than sparse matrix operations
	 * ------------------------------------
	 */
	RealType discount = 1;

	/*
  % ------------------------
  % method 1:
  %
  % BY(ib,ja) = B(ib,jb)*Y( jb,ja ),
  %
  % BY has size  nrow_B by ncol_A
  %
  % X(ib,ia) = BY(ib,ja)*At(ja,ia)
  % ------------------------
  */
	discount          = (is_dense_b) ? denseFlopDiscount : 1;
	RealType flops_by = discount * 2.0 * nnz_b * ncol_A;

	discount             = (is_dense_a) ? denseFlopDiscount : 1;
	RealType flops_by_at = discount * 2.0 * nnz_a * nrow_B;

	RealType flops_method1 = flops_by + flops_by_at;
	RealType nnz_method1   = nrow_B * ncol_A;

	/*
  % ------------------------
  % method 2:
  %
  % YAt(jb,ia)  = Y(jb,ja) * At(ja,ia)
  %
  % YAt has size ncol_B by nrow_A
  %
  % X(ib,ia) = B(ib,jb)*YAt(jb,ia)
  % ------------------------
  */
	discount            = (is_dense_a) ? denseFlopDiscount : 1;
	RealType flops_y_at = discount * 2.0 * nnz_a * ncol_B;

	discount    = (is_dense_b) ? denseFlopDiscount : 1;
	flops_by_at = discount * 2.0 * nnz_b * nrow_A;

	RealType flops_method2 = flops_y_at + flops_by_at;
	RealType nnz_method2   = ncol_B * nrow_A;

	/*
  % ------------------
  % method 3:
  % X = kron(A,B) * Y
  %
  % evaluate all nonzeros in A and nonzeros in B
  % need A and B to be very sparse
  % ------------------
  */
	RealType flops_method3 = 2.0 * nnz_a * nnz_b;
	RealType nnz_method3   = nnz_a * nnz_b;

	RealType kron_flops = std::min(flops_method1, std::min(flops_method2, flops_method3));

	RealType kron_nnz = std::min(nnz_method1, std::min(nnz_method2, nnz_method3));

	const bool minimize_flops = true;
	if (minimize_flops) {
		/*
		 * ----------------------------
		 * minimize the number of flops
		 * ----------------------------
		 */

		imethod = 3;
		if (kron_flops == flops_method1) {
			imethod = 1;
		};
		if (kron_flops == flops_method2) {
			imethod = 2;
		};
	} else {
		/*
		 * ------------------------------------
		 * minimize amount of temporary storage
		 * ------------------------------------
		 */
		imethod = 3;
		if (kron_nnz == nnz_method1) {
			imethod = 1;
		};
		if (kron_nnz == nnz_method2) {
			imethod = 2;
		};
	}

	*p_kron_nnz   = kron_nnz;
	*p_kron_flops = kron_flops;
	*p_imethod    = imethod;
}

#undef MIN
#undef TRUE
#undef FALSE
