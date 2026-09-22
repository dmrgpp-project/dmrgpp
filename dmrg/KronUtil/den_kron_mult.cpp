#include "util.h"

#include <Kokkos_Profiling_ScopedRegion.hpp>

template <typename ComplexOrRealType>
void denKronMultMethod(const int                                                   imethod,
                       const char                                                  transA,
                       const char                                                  transB,
                       const PsimagLite::Matrix<ComplexOrRealType>&                a_,
                       const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                       SizeType                                                    offsetY,
                       typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                       SizeType                                                    offsetX,
                       PsimagLite::GemmR<ComplexOrRealType>&                       gemmR)
{
	Kokkos::Profiling::ScopedRegion region("PsimagLite::csr_den_kron_mult_method");

	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int  nrow_a     = a_.n_row();
	const int  ncol_a     = a_.n_col();
	const int  nrow_b     = b_.n_row();
	const int  ncol_b     = b_.n_col();

	const int is_trans_a      = (transA == 'T') || (transA == 't');
	const int is_trans_b      = (transB == 'T') || (transB == 't');
	const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

	const int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	const int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;
	const int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	const int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

	const int nrow_x = nrow_2;
	const int ncol_x = nrow_1;
	const int nrow_y = ncol_2;
	const int ncol_y = ncol_1;

	PsimagLite::MatrixNonOwned<ComplexOrRealType>       xout(nrow_x, ncol_x, xout_, offsetX);
	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_y, ncol_y, yin_, offsetY);

	assert((imethod == 1) || (imethod == 2) || (imethod == 3));

	/*
	 *   -------------------------------------------------------------
	 *   A and B in dense matrix format
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

	if (imethod == 1) {
		Kokkos::Profiling::ScopedRegion region(
		    "PsimagLite::csr_den_kron_mult_method::imethod1");

		/*
		 *  --------------------------------------------
		 *  BY(iby,jby) = op(B(ib,jb))*Y(iy,jy)
		 *
		 *  X(ix,jx) += BY(iby,jby ) * transpose(op(A(ia,ja)))
		 *  --------------------------------------------
		 */
		const int                                           nrow_by = nrow_x;
		const int                                           ncol_by = ncol_y;
		PsimagLite::Matrix<ComplexOrRealType>               by(nrow_by, ncol_by);
		PsimagLite::MatrixNonOwned<ComplexOrRealType>       by_ref(by);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> by_const_ref(by);

		/*
		 * ---------------
		 * setup BY(ib,ja)
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
			 * BY(iby,jby)  = op(B(ib,jb))*Y(iy,jy)
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
			 * X(ix,jx) += BY(iby,jby) * transpose(op(A(ia,ja)))
			 * -------------------------------------------
			 */
			const char trans = (is_trans_a || is_conj_trans_a) ? 'N' : 'T';
			if (is_complex && is_conj_trans_a) {
				// --------------------------------------------
				// transpose( conj( transpose(A) ) ) is conj(A)
				// perform  conj operation
				// --------------------------------------------
				PsimagLite::Matrix<ComplexOrRealType> a_conj(nrow_a, ncol_a);

				for (int ja = 0; ja < ncol_a; ja++) {
					for (int ia = 0; ia < nrow_a; ia++) {
						a_conj(ia, ja) = PsimagLite::conj(a_(ia, ja));
					};
				};

				den_matmul_post(trans,
				                nrow_a,
				                ncol_a,
				                a_conj,

				                nrow_by,
				                ncol_by,
				                by_const_ref,

				                nrow_x,
				                ncol_x,
				                xout,
				                gemmR);
			} else {

				den_matmul_post(trans,
				                nrow_a,
				                ncol_a,
				                a_,

				                nrow_by,
				                ncol_by,
				                by_const_ref,

				                nrow_x,
				                ncol_x,
				                xout,
				                gemmR);
			};
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
		const int                                           nrow_y_at = nrow_y;
		const int                                           ncol_y_at = ncol_x;
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
			 * YAt(jb,ia) = Y(jb,ja) * tranpose(op(A(ia,ja)))
			 * ---------------------
			 */
			const char trans = (is_trans_a || is_conj_trans_a) ? 'N' : 'T';
			if (is_complex && is_conj_trans_a) {
				// --------------------------------------------
				// transpose( conj( transpose(A) ) ) is conj(A)
				// perform in-place conj operation
				// --------------------------------------------
				PsimagLite::Matrix<ComplexOrRealType> a_conj(nrow_a, ncol_a);

				for (int ja = 0; ja < ncol_a; ja++) {
					for (int ia = 0; ia < nrow_a; ia++) {
						a_conj(ia, ja) = PsimagLite::conj(a_(ia, ja));
					};
				};

				den_matmul_post(trans,
				                nrow_a,
				                ncol_a,
				                a_conj,

				                nrow_y,
				                ncol_y,
				                yin,

				                nrow_y_at,
				                ncol_y_at,
				                yat_ref,
				                gemmR);
			} else {

				den_matmul_post(trans,
				                nrow_a,
				                ncol_a,
				                a_,

				                nrow_y,
				                ncol_y,
				                yin,

				                nrow_y_at,
				                ncol_y_at,
				                yat_ref,
				                gemmR);
			};
		}

		{
			/*
			 * ------------
			 * X(ib,ia) += op(B(ib,jb)) * YAt(jb,ia)
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
		 *
		 * C = kron(transpose(A),B)
		 * C([ib,ja], [jb,ia]) = At(ja,ia) * B(ib,jb)
		 * X([ib,ja]) = B(ib,jb) * Y(jb,ia) * transpose(At(ja,ia))
		 * X([ib,ja]) = B(ib,jb) * Y(jb,ia) * A(ia,ja)
		 *            = (A(ia,ja)*B(ib,jb)) * Y(jb,ia)
		 *
		 * C = kron(A, transpose(B))
		 * C([jb,ia],[ib,ja]) = A(ia,ja) * Bt(jb,ib)
		 * X(jb,ia) = (A(ia,ja) * Bt(jb,ib)) * Y(ib,ja)
		 * X(jb,ia) = Bt(jb,ib) * Y(ib,ja) * transpose(A(ia,ja))
		 *          = Bt(jb,ib) * Y(ib,ja) * At(ja,ia)
		 *          = B(ib,jb)  * Y(ib,ja) * A(ia,ja)
		 *          = (A(ia,ja)*B(ib,jb)) * Y(ib,ja)
		 *
		 *
		 * C = kron( transpose(A), transpose(B))
		 * C([jb,ja], [ib,ia] ) = At(ja,ia) * Bt(jb,ib)
		 * X(jb,ja) = ( At(ja,ia) * Bt(jb,ib) ) * Y(ib,ia)
		 *          = Bt(jb,ib) * Y(ib,ia) * transpose(At(ja,ia))
		 *          = B(ib,jb) * Y(ib,ia) * A(ia,ja)
		 * ---------------------------------------------
		 */

		int ia = 0;
		int ja = 0;
		int ib = 0;
		int jb = 0;

		for (ia = 0; ia < nrow_a; ia++) {
			for (ja = 0; ja < ncol_a; ja++) {
				for (ib = 0; ib < nrow_b; ib++) {
					for (jb = 0; jb < ncol_b; jb++) {
						ComplexOrRealType aij = a_(ia, ja);
						if (is_complex && is_conj_trans_a) {
							aij = PsimagLite::conj(aij);
						};

						ComplexOrRealType bij = b_(ib, jb);
						if (is_complex && is_conj_trans_b) {
							bij = PsimagLite::conj(bij);
						};

						ComplexOrRealType cij = aij * bij;

						int ix = (is_trans_b || is_conj_trans_b) ? jb : ib;
						int jx = (is_trans_a || is_conj_trans_a) ? ja : ia;
						int iy = (is_trans_b || is_conj_trans_b) ? ib : jb;
						int jy = (is_trans_a || is_conj_trans_a) ? ia : ja;

						ComplexOrRealType yij = yin(iy, jy);
						xout(ix, jx) += (cij * yij);
					};
				};
			};
		};
	};
}

template <typename ComplexOrRealType>
void denKronMult(const char                                                  transA,
                 const char                                                  transB,
                 const PsimagLite::Matrix<ComplexOrRealType>&                a_,
                 const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                 SizeType                                                    offsetY,
                 typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                 SizeType                                                    offsetX,
                 const typename PsimagLite::Real<ComplexOrRealType>::Type    denseFlopDiscount,
                 PsimagLite::GemmR<ComplexOrRealType>&                       gemmR)
{
	/*
	 *   -------------------------------------------------------------
	 *   A and B in dense matrix format
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
	const int nrow_a = a_.n_row();
	const int ncol_a = a_.n_col();
	const int nrow_b = b_.n_row();
	const int ncol_b = b_.n_col();
	int       nnz_a  = nrow_a * ncol_a;
	int       nnz_b  = nrow_b * ncol_b;

	ComplexOrRealType kron_nnz   = 0;
	ComplexOrRealType kron_flops = 0;
	int               imethod    = 1;

	const int is_trans_a      = (transA == 'T') || (transA == 't');
	const int is_trans_b      = (transB == 'T') || (transB == 't');
	const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

	int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;

	int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

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

	den_kron_mult_method(imethod, transA, transB, a_, b_, yin, offsetY, xout, offsetX, gemmR);
}
#undef BY
#undef YAt
#undef X
#undef Y
#undef A
#undef B
#undef X2
