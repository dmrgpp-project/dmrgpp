#include "KronUtil.h"
#include "util.h"

#include <Kokkos_Core.hpp>
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

template <typename T> void runKronChecks()
{
	using RealT                     = typename PsimagLite::Real<T>::Type;
	const RealT dense_flop_discount = 0.2;

	static const bool    needs_printing     = false;
	const SizeType       gemm_rnb           = 49;
	const SizeType       threads_for_gemm_r = 1;
	PsimagLite::GemmR<T> gemm_r(needs_printing, gemm_rnb, threads_for_gemm_r);

	for (int threshold_b_idx = 0; threshold_b_idx <= 11; ++threshold_b_idx) {
		double threshold_b = .1 * threshold_b_idx;
		for (int threshold_a_idx = 0; threshold_a_idx <= 11; ++threshold_a_idx) {
			double threshold_a = .1 * threshold_a_idx;
			for (int ncol_a = 1; ncol_a <= 7; ncol_a += 3) {
				for (int nrow_a = 1; nrow_a <= 7; nrow_a += 3) {
					for (int ncol_b = 1; ncol_b <= 7; ncol_b += 3) {
						for (int nrow_b = 1; nrow_b <= 10; nrow_b += 3) {
							for (int itrans_a = 0; itrans_a <= 2;
							     itrans_a++) {
								for (int itrans_b = 0;
								     itrans_b <= 2;
								     itrans_b++) {
									char trans_a
									    = (itrans_a == 1)
									    ? 'T'
									    : ((itrans_a == 2)
									           ? 'C'
									           : 'N');
									char trans_b
									    = (itrans_b == 1)
									    ? 'T'
									    : ((itrans_b == 2)
									           ? 'C'
									           : 'N');

									int imethod = 0;

									int is_trans_a
									    = (trans_a == 'T');
									int is_trans_b
									    = (trans_b == 'T');
									int is_conj_trans_a
									    = (trans_a == 'C');
									int is_conj_trans_b
									    = (trans_b == 'C');

									int nrow_1
									    = (is_trans_a
									       || is_conj_trans_a)
									    ? ncol_a
									    : nrow_a;
									int ncol_1
									    = (is_trans_a
									       || is_conj_trans_a)
									    ? nrow_a
									    : ncol_a;

									int nrow_2
									    = (is_trans_b
									       || is_conj_trans_b)
									    ? ncol_b
									    : nrow_b;
									int ncol_2
									    = (is_trans_b
									       || is_conj_trans_b)
									    ? nrow_b
									    : ncol_b;

									int nrow_c
									    = nrow_1 * nrow_2;
									int ncol_c
									    = ncol_1 * ncol_2;

									int nrow_x = nrow_2;
									int ncol_x = nrow_1;

									int nrow_y = ncol_2;
									int ncol_y = ncol_1;

									PsimagLite::Matrix<T> a(
									    nrow_a, ncol_a);
									PsimagLite::Matrix<T> b(
									    nrow_b, ncol_b);

									PsimagLite::Matrix<T> y(
									    nrow_y, ncol_y);
									PsimagLite::MatrixNonOwned<
									    const T>
									    y_ref(y);

									PsimagLite::Matrix<T> x1(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    x1_ref(x1);
									PsimagLite::Matrix<T> x2(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    x2_ref(x2);
									PsimagLite::Matrix<T> x3(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    x3_ref(x3);
									PsimagLite::Matrix<T> x4(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    x4_ref(x4);

									PsimagLite::Matrix<T> sx1(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx1_ref(sx1);
									PsimagLite::Matrix<T> sx2(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx2_ref(sx2);
									PsimagLite::Matrix<T> sx3(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx3_ref(sx3);
									PsimagLite::Matrix<T> sx4(
									    nrow_x, ncol_x);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx4_ref(sx4);

									if (threshold_a == 0) {
										if (nrow_a
										    == ncol_a) {
											den_eye(
											    nrow_a,
											    ncol_a,
											    a);
											REQUIRE(
											    den_is_eye(
											        a));
											PsimagLite::
											    CrsMatrix<
											        T>
											        a(a);
											REQUIRE(
											    csr_is_eye(
											        a));
										} else {
											den_zeros(
											    nrow_a,
											    ncol_a,
											    a);
											REQUIRE(
											    den_is_zeros(
											        a));
											PsimagLite::
											    CrsMatrix<
											        T>
											        a(a);
											REQUIRE(
											    csr_is_zeros(
											        a));
										}
									} else {
										den_gen_matrix(
										    nrow_a,
										    ncol_a,
										    threshold_a,
										    a);
										PsimagLite::
										    CrsMatrix<T>
										        a(a);
										REQUIRE(
										    den_is_eye(a)
										    == csr_is_eye(
										        a));
										REQUIRE(
										    den_is_zeros(a)
										    == csr_is_zeros(
										        a));
									}

									if ((threshold_b == 0)
									    && (nrow_b == ncol_b)) {
										den_eye(nrow_b,
										        ncol_b,
										        b);
										REQUIRE(
										    den_is_eye(b));
										PsimagLite::
										    CrsMatrix<T>
										        b(b);
										REQUIRE(
										    csr_is_eye(b));
									} else {
										den_gen_matrix(
										    nrow_b,
										    ncol_b,
										    threshold_b,
										    b);
										PsimagLite::
										    CrsMatrix<T>
										        b(b);
										REQUIRE(
										    den_is_eye(b)
										    == csr_is_eye(
										        b));
										REQUIRE(
										    den_is_zeros(b)
										    == csr_is_zeros(
										        b));
									}

									den_gen_matrix(
									    nrow_y, ncol_y, 1.0, y);

									den_zeros(
									    nrow_x, ncol_x, x1);
									den_zeros(
									    nrow_x, ncol_x, x2);
									den_zeros(
									    nrow_x, ncol_x, x3);

									den_zeros(
									    nrow_x, ncol_x, sx1);
									den_zeros(
									    nrow_x, ncol_x, sx2);
									den_zeros(
									    nrow_x, ncol_x, sx3);

									imethod = 1;
									den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    x1_ref.getVector(),
									    0,
									    gemm_r);

									imethod = 2;
									den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    x2_ref.getVector(),
									    0,
									    gemm_r);
									imethod = 3;
									den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    x3_ref.getVector(),
									    0,
									    gemm_r);

									// ------------------
									// form C = kron(A,B)
									// ------------------
									PsimagLite::Matrix<T> c(
									    nrow_c, ncol_c);

									den_kron_form_general(
									    trans_a,
									    trans_b,
									    nrow_a,
									    ncol_a,
									    a,
									    nrow_b,
									    ncol_b,
									    b,
									    c);

									// -----------------------
									// perform matrix-multiply
									// -----------------------
									{
										const char trans1
										    = 'N';
										const char trans2
										    = 'N';
										const T alpha = 1.0;
										const T beta  = 0.0;

										// ------------------------------
										// reshape X, Y as
										// column vectors
										// ------------------------------
										const int mm
										    = nrow_x
										    * ncol_x;
										const int nn = 1;
										const int kk
										    = ncol_c;

										const int ld1
										    = nrow_c;
										const int ld2
										    = nrow_y
										    * ncol_y;
										const int ld3
										    = nrow_x
										    * ncol_x;

										const T* const p_a
										    = &(c(0, 0));
										const T* const p_b = &(
										    y_ref
										        .getVector()
										            [0]);
										T* p_c = &(
										    x4_ref
										        .getVector()
										            [0]);
										psimag::BLAS::GEMM(
										    trans1,
										    trans2,
										    mm,
										    nn,
										    kk,
										    alpha,
										    pA,
										    ld1,
										    pB,
										    ld2,
										    beta,
										    pC,
										    ld3);
									}

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff12
											    = std::abs(
											        x1(ix,
											           jx)
											        - x2(
											            ix,
											            jx));
											auto diff23
											    = std::abs(
											        x2(ix,
											           jx)
											        - x3(
											            ix,
											            jx));
											auto diff31
											    = std::abs(
											        x3(ix,
											           jx)
											        - x1(
											            ix,
											            jx));
											auto diff41
											    = std::abs(
											        x4(ix,
											           jx)
											        - x1(
											            ix,
											            jx));
											auto diffmax = std::max(
											    diff41,
											    std::max(
											        diff12,
											        std::max(
											            diff23,
											            diff31)));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den: transA="
												    << trans_a
												    << ", itransA "
												    << itrans_a
												    << " transB="
												    << trans_b
												    << ", itransB "
												    << itrans_b
												    << " nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff12 "
												    << diff12
												    << ", diff23 "
												    << diff23
												    << ", diff31 "
												    << diff31
												    << " diff41 "
												    << diff41
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									/*
									 * ------------------
									 * test sparse matrix
									 * ------------------
									 */
									PsimagLite::CrsMatrix<T> a(
									    a);
									REQUIRE(den_is_eye(a)
									        == csr_is_eye(a));
									REQUIRE(den_is_zeros(a)
									        == csr_is_zeros(a));
									PsimagLite::CrsMatrix<T> b(
									    b);
									REQUIRE(den_is_eye(b)
									        == csr_is_eye(b));
									REQUIRE(den_is_zeros(b)
									        == csr_is_zeros(b));

									imethod = 1;
									csr_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,

									    b,

									    y_ref,
									    sx1_ref);

									imethod = 2;
									csr_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,

									    b,

									    y_ref,
									    sx2_ref);

									imethod = 3;
									csr_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,

									    b,

									    y_ref,
									    sx3_ref);

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff1 = std::abs(
											    x1(ix,
											       jx)
											    - sx1(
											        ix,
											        jx));
											auto diff2 = std::abs(
											    x2(ix,
											       jx)
											    - sx2(
											        ix,
											        jx));
											auto diff3 = std::abs(
											    x3(ix,
											       jx)
											    - sx3(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "csr: transA="
												    << trans_a
												    << ", itransA "
												    << itrans_a
												    << " transB="
												    << trans_b
												    << ", itransB "
												    << itrans_b
												    << " nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff1 "
												    << diff1
												    << ", diff2 "
												    << diff2
												    << ", diff3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									/*
									 * ---------------------
									 * test generic interface
									 * ---------------------
									 */

									den_zeros(
									    nrow_x, ncol_x, x1);
									den_zeros(
									    nrow_x, ncol_x, sx1);

									den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    x1_ref.getVector(),
									    0,
									    gemm_r);

									csr_kron_mult(
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    sx1_ref.getVector(),
									    0,
									    RealT(
									        dense_flop_discount));

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff = std::abs(
											    x1(ix,
											       jx)
											    - sx1(
											        ix,
											        jx));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diff
											    > tol) {
												INFO(
												    "nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff "
												    << diff
												    << '\n');
												REQUIRE(
												    diff
												    <= tol);
											}
										}
									}

									/*
									 * -----------------------
									 * test mixed matrix types
									 * dense and CSR
									 * -----------------------
									 */

									den_zeros(
									    nrow_x, ncol_x, sx1);
									den_csr_kron_mult(

									    trans_a,
									    trans_b,

									    a,

									    b,

									    y_ref.getVector(),
									    0,
									    sx1_ref.getVector(),
									    0,
									    RealT(
									        dense_flop_discount),
									    gemm_r);

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff1 = std::abs(
											    x1(ix,
											       jx)
											    - sx1(
											        ix,
											        jx));
											auto diff2
											    = 0.0;
											auto diff3
											    = 0.0;
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: itransA "
												    << itrans_a
												    << "itransB "
												    << itrans_b
												    << " nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff1 "
												    << diff1
												    << ", diff2 "
												    << diff2
												    << ", diff3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									den_zeros(
									    nrow_x, ncol_x, sx1);
									den_zeros(
									    nrow_x, ncol_x, sx2);
									den_zeros(
									    nrow_x, ncol_x, sx3);

									imethod = 1;
									den_csr_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,

									    a,

									    b,

									    y_ref.getVector(),
									    0,
									    sx1_ref.getVector(),
									    0,
									    gemm_r);

									imethod = 2;
									den_csr_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,

									    a,

									    b,

									    y_ref.getVector(),
									    0,
									    sx2_ref.getVector(),
									    0,
									    gemm_r);

									imethod = 3;
									den_csr_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    sx3_ref.getVector(),
									    0,
									    gemm_r);

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff1 = std::abs(
											    x1(ix,
											       jx)
											    - sx1(
											        ix,
											        jx));
											auto diff2 = std::abs(
											    x2(ix,
											       jx)
											    - sx2(
											        ix,
											        jx));
											auto diff3 = std::abs(
											    x3(ix,
											       jx)
											    - sx3(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: itransA "
												    << itrans_a
												    << "itransB "
												    << itrans_b
												    << " nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff1 "
												    << diff1
												    << ", diff2 "
												    << diff2
												    << ", diff3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									/*
									 * -----------------------
									 * test mixed matrix types
									 * CSR and dense
									 * -----------------------
									 */
									den_zeros(
									    nrow_x, ncol_x, sx1);

									csr_den_kron_mult(
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    0,
									    sx1_ref.getVector(),
									    0,
									    RealT(
									        dense_flop_discount),
									    gemm_r);

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff1 = std::abs(
											    x1(ix,
											       jx)
											    - sx1(
											        ix,
											        jx));
											auto diff2
											    = 0.0;
											auto diff3
											    = 0.0;
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: itransA "
												    << itrans_a
												    << "itransB "
												    << itrans_b
												    << " nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff1 "
												    << diff1
												    << ", diff2 "
												    << diff2
												    << ", diff3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									den_zeros(
									    nrow_x, ncol_x, sx1);
									den_zeros(
									    nrow_x, ncol_x, sx2);
									den_zeros(
									    nrow_x, ncol_x, sx3);

									imethod              = 1;
									const SizeType izero = 0;
									csr_den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    izero,
									    sx1_ref.getVector(),
									    izero,
									    gemm_r);

									imethod = 2;
									csr_den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    izero,
									    sx2_ref.getVector(),
									    izero,
									    gemm_r);

									imethod = 3;
									csr_den_kron_mult_method(
									    imethod,
									    trans_a,
									    trans_b,
									    a,
									    b,
									    y_ref.getVector(),
									    izero,
									    sx3_ref.getVector(),
									    izero,
									    gemm_r);

									for (int jx = 0;
									     jx < ncol_x;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_x;
										     ++ix) {
											auto diff1 = std::abs(
											    x1(ix,
											       jx)
											    - sx1(
											        ix,
											        jx));
											auto diff2 = std::abs(
											    x2(ix,
											       jx)
											    - sx2(
											        ix,
											        jx));
											auto diff3 = std::abs(
											    x3(ix,
											       jx)
											    - sx3(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: itransA "
												    << itrans_a
												    << "itransB "
												    << itrans_b
												    << " nrow_A "
												    << nrow_a
												    << " ncol_A "
												    << ncol_a
												    << " nrow_B "
												    << nrow_b
												    << " ncol_B "
												    << ncol_b
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff1 "
												    << diff1
												    << ", diff2 "
												    << diff2
												    << ", diff3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

TEST_CASE("kron_mult_test1_double", "[kron][basic]") { runKronChecks<double>(); }

TEST_CASE("kron_mult_test1_complex", "[kron][basic]") { runKronChecks<std::complex<double>>(); }

int main(int argc, char* argv[])
{
	Kokkos::initialize(argc, argv);
	int result = Catch::Session().run(argc, argv);
	Kokkos::finalize();
	return result;
}
