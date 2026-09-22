#include "util.h"

#include <Kokkos_Profiling_ScopedRegion.hpp>

template <typename ComplexOrRealType>
void csrDenKronMultMethod(const int                                                   imethod,
                          const char                                                  transA,
                          const char                                                  transB,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&             a_,
                          const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                          SizeType                                                    offsetY,
                          typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                          SizeType                                                    offsetX,
                          PsimagLite::GemmR<ComplexOrRealType>&                       gemmR)
{
	Kokkos::Profiling::ScopedRegion region("PsimagLite::csr_den_kron_mult_method");

	const bool is_complex      = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int  is_trans_a      = (transA == 'T') || (transA == 't');
	const int  is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int  is_trans_b      = (transB == 'T') || (transB == 't');
	const int  is_conj_trans_b = (transB == 'C') || (transB == 'c');

	const int nrow_a = a_.rows();
	const int ncol_a = a_.cols();
	const int nrow_b = b_.rows();
	const int ncol_b = b_.cols();

	const int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	const int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;
	const int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	const int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

	const int nrow_x = nrow_2;
	const int ncol_x = nrow_1;
	const int nrow_y = ncol_2;
	const int ncol_y = ncol_1;

	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_y, ncol_y, yin_, offsetY);
	PsimagLite::MatrixNonOwned<ComplexOrRealType>       xout(nrow_x, ncol_x, xout_, offsetX);

	assert((imethod == 1) || (imethod == 2) || (imethod == 3));

	bool no_work = (csr_is_zeros(a_) || den_is_zeros(b_));
	if (no_work) {
		return;
	};
	/*
	 *   -------------------------------------------------------------
	 *   A in compressed sparse ROW format
	 *   B in dense format
	 *
	 *   X += kron( op(A), op(B)) * Y
	 *   X +=  op(B) * Y * transpose(op(A))
	 *
	 *   nrow_X = nrow_2,   ncol_X = nrow_1
	 *   nrow_Y = ncol_2,   nrow_Y = nrow_2
	 *
	 *   that can be computed as either
	 *   imethod == 1
	 *
	 *   X(ib,ia) +=  (B(ib,jb) * Y( jb,ja)) * transpose( A(ia,ja))
	 *
	 *   X(ix,jx) +=  (B2(ib2,jb2) * Y(iy,jy)) * transpose(A1(ia1,ja1))
	 *
	 *
	 *
	 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
	 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
	 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops
	 *
	 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
	 *
	 *   imethod == 2
	 *
	 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
	 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))
	 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
	 *
	 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
	 *
	 *   imethod == 3
	 *
	 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B
	 *
	 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
	 *   -------------------------------------------------------------
	 */

	if (imethod == 1) {
		Kokkos::Profiling::ScopedRegion region(
		    "PsimagLite::csr_den_kron_mult_method::imethod1");

		/*
		 *  --------------------------------------------
		 *  BY(ib,ja) = (B(ib,jb))*Y(jb,ja)
		 *
		 *  X(ib,ia) += BY(ib,ja ) * transpose(A(ia,ja))
		 *  --------------------------------------------
		 */

		int                                                 nrow_by = nrow_x;
		int                                                 ncol_by = ncol_y;
		PsimagLite::Matrix<ComplexOrRealType>               by(nrow_by, ncol_by);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> by_const_ref(by);
		PsimagLite::MatrixNonOwned<ComplexOrRealType>       by_ref(by);

		/*
		 * ---------------
		 * setup BY
		 * ---------------
		 */

		{
			int iby = 0;
			int jby = 0;

			// not needed FIXME
			for (jby = 0; jby < ncol_by; jby++) {
				for (iby = 0; iby < nrow_by; iby++) {
					by(iby, jby) = 0;
				};
			};
		}

		{
			/*
			 * ------------------------------
			 * BY(ib,ja)  = B(ib,jb)*Y(jb,ja)
			 * ------------------------------
			 */
			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;

			den_matmul_pre(trans,
			               nrow_b,
			               ncol_b,
			               b_,

			               nrow_y,
			               ncol_y,
			               yin,

			               nrow_by,
			               ncol_by,
			               by_ref,
			               gemmR);
		}

		{
			/*
			 * -------------------------------------------
			 * X(ib,ia) += BY(ib,ja) * transpose(A(ia,ja))
			 * -------------------------------------------
			 */

			/*
			 * ----------------------------------
			 * Note trans = 'Z' means use conj(A)
			 * ----------------------------------
			 */
			const char trans = is_trans_a ? 'N' : (is_conj_trans_a ? 'Z' : 'T');
			csr_matmul_post(trans,
			                a_,

			                nrow_by,
			                ncol_by,
			                by_const_ref,

			                nrow_x,
			                ncol_x,
			                xout);
		}
	} else if (imethod == 2) {
		Kokkos::Profiling::ScopedRegion region(
		    "PsimagLite::csr_den_kron_mult_method::imethod2");

		/*
		 * ---------------------
		 * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja))
		 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
		 * ---------------------
		 */

		int                                                 nrow_y_at = nrow_y;
		int                                                 ncol_y_at = ncol_x;
		PsimagLite::Matrix<ComplexOrRealType>               yat(nrow_y_at, ncol_y_at);
		PsimagLite::MatrixNonOwned<ComplexOrRealType>       yat_ref(yat);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> yat_const_ref(yat);

		/*
		 * ----------------
		 * setup YAt(jb,ia)
		 * ----------------
		 */

		{
			int iy = 0;
			int jy = 0;

			// not needed, FIXME
			for (jy = 0; jy < ncol_y_at; jy++) {
				for (iy = 0; iy < nrow_y_at; iy++) {
					yat(iy, jy) = 0;
				};
			};
		}

		{
			/*
			 * ---------------------
			 * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja)
			 * ---------------------
			 */
			/*
			 * ------------------------------
			 * note transa = 'Z' mean conj(A)
			 * ------------------------------
			 */
			const char transa = is_trans_a ? 'N' : (is_conj_trans_a ? 'Z' : 'T');
			csr_matmul_post(transa,
			                a_,

			                nrow_y,
			                ncol_y,
			                yin,

			                nrow_y_at,
			                ncol_y_at,
			                yat_ref);
		}

		{
			/*
			 * ------------
			 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
			 * ------------
			 */

			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;
			den_matmul_pre(trans,
			               nrow_b,
			               ncol_b,
			               b_,

			               nrow_y_at,
			               ncol_y_at,
			               yat_const_ref,

			               nrow_x,
			               ncol_x,
			               xout,
			               gemmR);
		}
	} else if (imethod == 3) {
		Kokkos::Profiling::ScopedRegion region(
		    "PsimagLite::csr_den_kron_mult_method::imethod3");

		/*
		 * ---------------------------------------------
		 * C = kron(A,B)
		 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
		 * X([ib,ia]) += C([ib,ia],[jb,ja]) * Y([jb,ja])
		 * ---------------------------------------------
		 */
		const bool use_inner_sparse = false;

		if (use_inner_sparse) {

			int ib = 0;
			int jb = 0;

			for (ib = 0; ib < nrow_b; ib++) {
				for (jb = 0; jb < ncol_b; jb++) {
					ComplexOrRealType bij = b_(ib, jb);
					if (is_complex && is_conj_trans_b) {
						bij = PsimagLite::conj(bij);
					};

					int ia = 0;
					for (ia = 0; ia < nrow_a; ia++) {
						int istarta = a_.getRowPtr(ia);
						int ienda   = a_.getRowPtr(ia + 1);

						int ka = 0;
						for (ka = istarta; ka < ienda; ka++) {
							int               ja  = a_.getCol(ka);
							ComplexOrRealType aij = a_.getValue(ka);
							if (is_complex && is_conj_trans_a) {
								aij = PsimagLite::conj(aij);
							};

							ComplexOrRealType cij = aij * bij;

							int ix = (is_trans_b || is_conj_trans_b)
							    ? jb
							    : ib;
							int jx = (is_trans_a || is_conj_trans_a)
							    ? ja
							    : ia;
							int iy = (is_trans_b || is_conj_trans_b)
							    ? ib
							    : jb;
							int jy = (is_trans_a || is_conj_trans_a)
							    ? ia
							    : ja;

							xout(ix, jx) += (cij * yin(iy, jy));
						};
					};
				};
			};
		} else {

			int ia = 0;
			for (ia = 0; ia < nrow_a; ia++) {
				int istarta = a_.getRowPtr(ia);
				int ienda   = a_.getRowPtr(ia + 1);
				int ka      = 0;
				for (ka = istarta; ka < ienda; ka++) {
					int               ja  = a_.getCol(ka);
					ComplexOrRealType aij = a_.getValue(ka);
					if (is_complex && is_conj_trans_a) {
						aij = PsimagLite::conj(aij);
					};

					int ib = 0;
					int jb = 0;

					for (ib = 0; ib < nrow_b; ib++) {
						for (jb = 0; jb < ncol_b; jb++) {
							ComplexOrRealType bij = b_(ib, jb);
							if (is_complex && is_conj_trans_b) {
								bij = PsimagLite::conj(bij);
							};

							ComplexOrRealType cij = aij * bij;

							int ix = (is_trans_b || is_conj_trans_b)
							    ? jb
							    : ib;
							int jx = (is_trans_a || is_conj_trans_a)
							    ? ja
							    : ia;
							int iy = (is_trans_b || is_conj_trans_b)
							    ? ib
							    : jb;
							int jy = (is_trans_a || is_conj_trans_a)
							    ? ia
							    : ja;

							xout(ix, jx) += (cij * yin(iy, jy));
						};
					};
				};
			};
		};
	};
}

template <typename ComplexOrRealType>
void csrDenKronMult(const char                                                  transA,
                    const char                                                  transB,
                    const PsimagLite::CrsMatrix<ComplexOrRealType>&             a_,
                    const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                    SizeType                                                    offsetY,
                    typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                    SizeType                                                    offsetX,
                    const typename PsimagLite::Real<ComplexOrRealType>::Type    denseFlopDiscount,
                    PsimagLite::GemmR<ComplexOrRealType>&                       gemmR)

{
	const int idebug = 0;
	/*
	 *   -------------------------------------------------------------
	 *   A in compressed sparse ROW format
	 *   B in dense matrix format
	 *
	 *   X += kron( A, B) * Y
	 *   that can be computed as either
	 *   imethod == 1
	 *
	 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
	 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
	 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops
	 *
	 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
	 *
	 *   imethod == 2
	 *
	 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
	 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))
	 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
	 *
	 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
	 *
	 *   imethod == 3
	 *
	 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B
	 *
	 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
	 *   -------------------------------------------------------------
	 */

	const int nrow_a = a_.rows();
	const int ncol_a = a_.cols();
	const int nrow_b = b_.rows();
	const int ncol_b = b_.cols();
	int       nnz_a  = csr_nnz(a_);
	int       nnz_b  = den_nnz(b_);

	bool no_work = (csr_is_zeros(a_) || den_is_zeros(b_));
	if (no_work) {
		if (idebug >= 1) {
			printf("csr_den: no_work, nrow_A=%d,ncol_A=%d,nrow_B=%d,ncol_B=%d\n",
			       nrow_a,
			       ncol_a,
			       nrow_b,
			       ncol_b);
		};

		return;
	};

	ComplexOrRealType kron_nnz   = 0;
	ComplexOrRealType kron_flops = 0;
	int               imethod    = 1;

	const int is_trans_a      = (transA == 'T') || (transA == 't');
	const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int is_trans_b      = (transB == 'T') || (transB == 't');
	const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

	/*
	 * -------------------------------------
	 * check for special case where A is eye
	 * -------------------------------------
	 */

	if (csr_is_eye(a_)) {
		if (idebug >= 1) {
			printf("csr_den: A is eye: nrow_A=%d,ncol_A=%d,nrow_B=%d,ncol_B=%d\n",
			       nrow_a,
			       ncol_a,
			       nrow_b,
			       ncol_b);
		};

		/*
	----------------------
	X +=  ( op(B) ) Y * transpose( op(A) )
	----------------------
	*/
		const int nrow_y = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;
		const int ncol_y = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;

		const int nrow_x = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
		const int ncol_x = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;

		PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(
		    nrow_y, ncol_y, yin_, offsetY);
		PsimagLite::MatrixNonOwned<ComplexOrRealType> xout(nrow_x, ncol_x, xout_, offsetX);

		// const char  trans1 =  (isTransB) ? 'T' : 'N';
		const char trans1 = transB;
		den_matmul_pre(
		    trans1, nrow_b, ncol_b, b_, nrow_y, ncol_y, yin, nrow_x, ncol_x, xout, gemmR);

		return;
	};

	const int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	const int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;

	const int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	const int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

	estimate_kron_cost(nrow_1,
	                   ncol_1,
	                   nnz_a,
	                   nrow_2,
	                   ncol_2,
	                   nnz_b,
	                   &kron_nnz,
	                   &kron_flops,
	                   &imethod,
	                   denseFlopDiscount);

	csr_den_kron_mult_method(
	    imethod, transA, transB, a_, b_, yin_, offsetY, xout_, offsetX, gemmR);
}

#undef B
#undef BY
#undef YAt
#undef X
#undef Y
