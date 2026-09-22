#include "util.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <PsimagLite/KokkosType.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

template <typename ComplexOrRealType>
void csrToDen(const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
              PsimagLite::Matrix<ComplexOrRealType>&          a_)
{
	int ia = 0;
	int ja = 0;

	int nrow_a = a.row();
	int ncol_a = a.col();

	for (ja = 0; ja < ncol_a; ja++) {
		for (ia = 0; ia < nrow_a; ia++) {
			a_(ia, ja) = 0;
		};
	};

	for (ia = 0; ia < nrow_a; ia++) {
		int istarta = a.getRowPtr(ia);
		int ienda   = a.getRowPtr(ia + 1);
		int ka      = 0;
		for (ka = istarta; ka < ienda; ka++) {
			ComplexOrRealType aij = a.getValue(ka);
			int               ja  = a.getCol(ka);
			a_(ia, ja)            = aij;
		};
	};
}

template <typename ComplexOrRealType>
void csrKronMultMethod(const int  imethod,
                       const char transA,
                       const char transB,

                       const PsimagLite::CrsMatrix<ComplexOrRealType>& a,

                       const PsimagLite::CrsMatrix<ComplexOrRealType>& b,

                       const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                       PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout)
{
	Kokkos::Profiling::ScopedRegion region("PsimagLite::csr_kron_mult_method");

	const int is_trans_a      = (transA == 'T') || (transA == 't');
	const int is_trans_b      = (transB == 'T') || (transB == 't');
	const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

	const int nrow_a = a.rows();
	const int ncol_a = a.cols();
	const int nrow_b = b.rows();
	const int ncol_b = b.cols();

	const int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	const int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;
	const int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	const int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

	const int nrow_x = nrow_2;
	const int ncol_x = nrow_1;
	const int nrow_y = ncol_2;
	const int ncol_y = ncol_1;

	assert((imethod == 1) || (imethod == 2) || (imethod == 3));

	bool no_work = (csr_is_zeros(a) || csr_is_zeros(b));
	if (no_work) {
		return;
	};
	/*
	 *   -------------------------------------------------------------
	 *   A and B in compressed sparse ROW format
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
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod1");

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
		PsimagLite::MatrixNonOwned<ComplexOrRealType>       by_ref(by);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> by_const_ref(by);
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
			csr_matmul_pre(trans,
			               b,

			               nrow_y,
			               ncol_y,
			               yin,

			               nrow_by,
			               ncol_by,
			               by_ref);
		}

		{
			/*
			 * -------------------------------------------
			 * X(ib,ia) += BY(ib,ja) * transpose(A(ia,ja))
			 * -------------------------------------------
			 */
			/*
			 * ---------------------------------
			 * note trans = 'Z' mean use conj(A)
			 * ---------------------------------
			 */
			const char trans = is_trans_a ? 'N' : (is_conj_trans_a ? 'Z' : 'T');
			csr_matmul_post(trans,
			                a,

			                nrow_by,
			                ncol_by,
			                by_const_ref,

			                nrow_x,
			                ncol_x,
			                xout);
		}
	} else if (imethod == 2) {
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod2");

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

			// not needed FIXME
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
			 * ---------------------------------
			 * note trans = 'Z' mean use conj(A)
			 * ---------------------------------
			 */
			const char transa = is_trans_a ? 'N' : (is_conj_trans_a ? 'Z' : 'T');
			csr_matmul_post(transa,
			                a,

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

			csr_matmul_pre(trans,
			               b,

			               nrow_y_at,
			               ncol_y_at,
			               yat_const_ref,

			               nrow_x,
			               ncol_x,
			               xout);
		}
	} else if (imethod == 3) {
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod3");

		/*
		 * ---------------------------------------------
		 * C = kron(A,B)
		 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
		 * X([ib,ia]) += C([ib,ia],[jb,ja]) * Y([jb,ja])
		 * ---------------------------------------------
		 */

		using ExecutionSpace = Kokkos::DefaultExecutionSpace;
		using MemorySpace    = ExecutionSpace::memory_space;
		using KokkosScalar   = typename PsimagLite::KokkosType<ComplexOrRealType>::type;

		int nnz_a = a.nonZeros();
		int nnz_b = b.nonZeros();

		Kokkos::View<int*, Kokkos::HostSpace> a_row_h(
		    Kokkos::view_alloc(Kokkos::WithoutInitializing,
		                       "PsimgLite:csr_kron_mult::imethod3::A_row_h"),
		    nnz_a);
		Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> a_col_h(
		    &a.getCol(0), nnz_a);
		Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    a_val_h(reinterpret_cast<const KokkosScalar*>(&a.getValue(0)), nnz_a);
		Kokkos::View<int*, Kokkos::HostSpace> b_row_h(
		    Kokkos::view_alloc(Kokkos::WithoutInitializing,
		                       "PsimgLite:csr_kron_mult::imethod3::B_row_h"),
		    nnz_b);
		Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> b_col_h(
		    &b.getCol(0), nnz_b);
		Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    b_val_h(reinterpret_cast<const KokkosScalar*>(&b.getValue(0)), nnz_b);

		{
			Kokkos::Profiling::ScopedRegion region(
			    "PsimgLite::csr_kron_mult_method::imethod3::store_rows");

			// store row for every entry
			int idx = 0;
			for (int ia = 0; ia < nrow_a; ++ia) {
				int istart = a.getRowPtr(ia);
				int iend   = a.getRowPtr(ia + 1);
				for (int ka = istart; ka < iend; ++ka) {
					a_row_h[idx] = ia;
					++idx;
				}
			}

			idx = 0;
			for (int ib = 0; ib < nrow_b; ++ib) {
				int istart = b.getRowPtr(ib);
				int iend   = b.getRowPtr(ib + 1);
				for (int kb = istart; kb < iend; ++kb) {
					b_row_h[idx] = ib;
					++idx;
				}
			}
		}

		auto a_row_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), a_row_h);
		auto a_col_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), a_col_h);
		auto a_val_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), a_val_h);
		auto b_row_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), b_row_h);
		auto b_col_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), b_col_h);
		auto b_val_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), b_val_h);

		auto yin_host = Kokkos::View<const KokkosScalar**,
		                             Kokkos::LayoutLeft,
		                             Kokkos::HostSpace,
		                             Kokkos::MemoryUnmanaged>(
		    reinterpret_cast<const KokkosScalar*>(&yin(0, 0)), nrow_y, ncol_y);

		auto y_dev = Kokkos::create_mirror_view_and_copy(
		    Kokkos::view_alloc(ExecutionSpace {}, MemorySpace {}), yin_host);

		auto x_dev = Kokkos::View<KokkosScalar**>(
		    "PsimgLite:csr_kron_mult::imethod3::x_dev", nrow_x, ncol_x);

		Kokkos::parallel_for(
		    "PsimgLite:csr_kron_mult::imethod3::spmv_kernel",
		    Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<2>>({ 0, 0 },
		                                                           { nnz_b, nnz_a }),
		    KOKKOS_LAMBDA(const size_t ib_idx, const size_t ia_idx) {
			    constexpr bool is_complex
			        = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
			    const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
			    const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

			    int          ia  = a_row_dev(ia_idx);
			    int          ja  = a_col_dev(ia_idx);
			    KokkosScalar aij = a_val_dev(ia_idx);
			    if constexpr (is_complex)
				    if (is_conj_trans_a)
					    aij = Kokkos::conj(aij);

			    int          ib  = b_row_dev(ib_idx);
			    int          jb  = b_col_dev(ib_idx);
			    KokkosScalar bij = b_val_dev(ib_idx);
			    if constexpr (is_complex)
				    if (is_conj_trans_b)
					    bij = Kokkos::conj(bij);

			    KokkosScalar cij = aij * bij;

			    int ix = (is_trans_b || is_conj_trans_b) ? jb : ib;
			    int jx = (is_trans_a || is_conj_trans_a) ? ja : ia;
			    int iy = (is_trans_b || is_conj_trans_b) ? ib : jb;
			    int jy = (is_trans_a || is_conj_trans_a) ? ia : ja;

			    KokkosScalar prod = cij * y_dev(iy, jy);
			    Kokkos::atomic_add(&x_dev(ix, jx), prod);
		    });
		{
			Kokkos::Profiling::ScopedRegion region(
			    "PsimgLite::csr_kron_mult_method::imethod3::copy_results");
			auto xhost
			    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace {}, x_dev);
			for (int ix = 0; ix < nrow_x; ++ix) {
				for (int jx = 0; jx < ncol_x; ++jx)
					xout(ix, jx)
					    += static_cast<ComplexOrRealType>(xhost(ix, jx));
			}
		}
	};
}

template <typename ComplexOrRealType>
void csrKronMultMethod(const int                                                   imethod,
                       const char                                                  transA,
                       const char                                                  transB,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>&             a,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>&             b,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                       SizeType                                                    offsetY,
                       typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                       SizeType                                                    offsetX)

{
	const int is_trans_a      = (transA == 'T') || (transA == 't');
	const int is_trans_b      = (transB == 'T') || (transB == 't');
	const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

	const int nrow_a = a.rows();
	const int ncol_a = a.cols();
	const int nrow_b = b.rows();
	const int ncol_b = b.cols();

	const int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	const int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;
	const int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	const int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

	const int                                           nrow_x = nrow_2;
	const int                                           ncol_x = nrow_1;
	const int                                           nrow_y = ncol_2;
	const int                                           ncol_y = ncol_1;
	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_y, ncol_y, yin_, offsetY);
	PsimagLite::MatrixNonOwned<ComplexOrRealType>       xout(nrow_x, ncol_x, xout_, offsetX);
	csr_kron_mult_method(imethod, transA, transB, a, b, yin, xout);
}

template <typename ComplexOrRealType>
void csrKronMult(const char                                                  transA,
                 const char                                                  transB,
                 const PsimagLite::CrsMatrix<ComplexOrRealType>&             a,
                 const PsimagLite::CrsMatrix<ComplexOrRealType>&             b,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                 SizeType                                                    offsetY,
                 typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                 SizeType                                                    offsetX,
                 const typename PsimagLite::Real<ComplexOrRealType>::Type    denseFlopDiscount)
{
	/*
	 *   -------------------------------------------------------------
	 *   A and B in compressed sparse ROW format
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
	int nnz_a = csrNnz(a);
	int nnz_b = csrNnz(b);

	bool no_work = (csrIsZeros(a) || csrIsZeros(b));
	if (no_work) {
		return;
	};

	ComplexOrRealType kron_nnz   = 0;
	ComplexOrRealType kron_flops = 0;
	int               imethod    = 1;

	const int is_trans_a      = (transA == 'T') || (transA == 't');
	const int is_trans_b      = (transB == 'T') || (transB == 't');
	const int is_conj_trans_a = (transA == 'C') || (transA == 'c');
	const int is_conj_trans_b = (transB == 'C') || (transB == 'c');

	const int nrow_a = a.rows();
	const int ncol_a = a.cols();
	const int nrow_b = b.rows();
	const int ncol_b = b.cols();

	// -----------------------------------
	// both A and B are considered sparse
	// -----------------------------------

	const int nrow_1 = (is_trans_a || is_conj_trans_a) ? ncol_a : nrow_a;
	const int ncol_1 = (is_trans_a || is_conj_trans_a) ? nrow_a : ncol_a;

	const int nrow_2 = (is_trans_b || is_conj_trans_b) ? ncol_b : nrow_b;
	const int ncol_2 = (is_trans_b || is_conj_trans_b) ? nrow_b : ncol_b;

	estimateKronCost(nrow_1,
	                 ncol_1,
	                 nnz_a,
	                 nrow_2,
	                 ncol_2,
	                 nnz_b,
	                 &kron_nnz,
	                 &kron_flops,
	                 &imethod,
	                 denseFlopDiscount);

	csrKronMultMethod(imethod, transA, transB, a, b, yin, offsetY, xout, offsetX);
}
