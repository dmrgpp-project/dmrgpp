#include "KronUtil.h"
#include "util.h"

#include <Kokkos_Core.hpp>
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

template <typename ComplexOrRealType> void runKronSubmatrixChecks()
{
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;

	for (int threshold_b_idx = 0; threshold_b_idx <= 11; ++threshold_b_idx) {
		double threshold_b = .1 * threshold_b_idx;
		for (int threshold_a_idx = 0; threshold_a_idx <= 11; ++threshold_a_idx) {
			double threshold_a = .1 * threshold_a_idx;
			for (int ncol_a = 1; ncol_a <= 10; ncol_a += 3) {
				for (int nrow_a = 1; nrow_a <= 10; nrow_a += 3) {
					for (int ncol_b = 1; ncol_b <= 10; ncol_b += 3) {
						for (int nrow_b = 1; nrow_b <= 10; nrow_b += 3) {
							PsimagLite::Matrix<ComplexOrRealType> a(
							    nrow_a, ncol_a);
							PsimagLite::Matrix<ComplexOrRealType> b(
							    nrow_b, ncol_b);

							if ((threshold_a == 0)
							    && (nrow_a == ncol_a)) {
								/*
								 * ------------------------------------
								 * special case to test identity
								 * matrix
								 * ------------------------------------
								 */
								den_eye(nrow_a, ncol_a, a);
								REQUIRE(den_is_eye(a));
								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    a(a);
								REQUIRE(csr_is_eye(a));
							} else {
								den_gen_matrix(
								    nrow_a, ncol_a, threshold_a, a);

								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    a(a);
								REQUIRE(den_is_eye(a)
								        == csr_is_eye(a));
								REQUIRE(den_is_zeros(a)
								        == csr_is_zeros(a));
							}

							if ((threshold_b == 0)
							    && (nrow_b == ncol_b)) {
								den_eye(nrow_b, ncol_b, b);
								REQUIRE(den_is_eye(b));
								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    b(b);
								REQUIRE(csr_is_eye(b));
							} else {
								den_gen_matrix(
								    nrow_b, ncol_b, threshold_b, b);

								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    b(b);
								REQUIRE(den_is_eye(b)
								        == csr_is_eye(b));
								REQUIRE(den_is_zeros(b)
								        == csr_is_zeros(b));
							}
							/*
							 * -----------------------------------
							 * generate compressed row format from
							 * dense matrices A, B
							 * -----------------------------------
							 */
							PsimagLite::CrsMatrix<ComplexOrRealType> a(
							    a);
							REQUIRE(den_is_eye(a) == csr_is_eye(a));

							PsimagLite::CrsMatrix<ComplexOrRealType> b(
							    b);
							REQUIRE(den_is_eye(a) == csr_is_eye(a));

							/*
							 * -----------------------------
							 * explicitly form C = kron(A,B)
							 * -----------------------------
							 */

							const int nrow_c = nrow_a * nrow_b;
							const int ncol_c = ncol_a * ncol_b;
							if (ncol_c < 2)
								continue;
							PsimagLite::Matrix<ComplexOrRealType> c(
							    nrow_c, ncol_c);

							den_kron_form(nrow_a,
							              ncol_a,
							              a,
							              nrow_b,
							              ncol_b,
							              b,
							              c);

							/*
							 * ---------------------------------------
							 * generate compressed sparse version of C
							 * ---------------------------------------
							 */
							PsimagLite::CrsMatrix<ComplexOrRealType> c(
							    c);

							PsimagLite::Vector<int>::Type rindex(
							    nrow_c);
							PsimagLite::Vector<int>::Type cindex(
							    ncol_c);

							/*
							 * --------------------
							 * extract  even rows
							 * extract odd columns
							 * --------------------
							 */
							int nrindex = 0;
							for (int ic = 0; ic < nrow_c; ic += 2) {
								rindex[nrindex++] = ic;
							}
							int ncindex = 0;
							for (int jc = 1; jc < ncol_c; jc += 2) {
								cindex[ncindex++] = jc;
							}

							/*
							 * -------------------------------
							 * extract submatrix from C into D
							 * -------------------------------
							 */

							int nrow_d = nrindex;
							int ncol_d = ncindex;
							PsimagLite::Matrix<ComplexOrRealType> d(
							    nrow_d, ncol_d);

							den_submatrix(nrow_c,
							              ncol_c,
							              c,
							              nrindex,
							              ncindex,
							              rindex,
							              cindex,
							              d);
							/*
							 * -------------------------------------
							 * generate submatrix from sparse version C
							 * -------------------------------------
							 */

							const int max_nnz_d = 1 + den_nnz(d);

							PsimagLite::CrsMatrix<ComplexOrRealType> sd(
							    d.n_row(), d.n_col());

							csr_submatrix(c,

							              nrindex,
							              ncindex,
							              max_nnz_d,
							              rindex,
							              cindex,

							              sd);

							/*
							 * -----------------------
							 * convert to dense matrix
							 * -----------------------
							 */
							PsimagLite::Matrix<ComplexOrRealType> dd(
							    nrow_d, ncol_d);
							crsMatrixToFullMatrix(dd, sd);

							/*
							 * --------------------------------------
							 * check both matrices should be the same
							 * --------------------------------------
							 */

							for (int jd = 0; jd < ncol_d; ++jd) {
								for (int id = 0; id < nrow_d;
								     ++id) {
									RealType diff = std::abs(
									    dd(id, jd) - d(id, jd));
									const RealType tol = 1.0
									    / (1000.0 * 1000.0);
									if (diff > tol) {
										INFO("nrow_D "
										     << nrow_d
										     << " ncol_D "
										     << ncol_d
										     << " DD(" << id
										     << "," << jd
										     << ")"
										     << dd(id, jd)
										     << " D(" << id
										     << "," << jd
										     << ")"
										     << d(id, jd)
										     << " di"
										        "ff "
										     << diff);
										REQUIRE(diff
										        <= tol);
									}
								}
							}

							/*
							 * --------------------------------
							 * form E = C(rindex(:),cindex(:))
							 * without forming C
							 * E should be the same as matrix D
							 * --------------------------------
							 */
							int nrow_e = nrindex;
							int ncol_e = ncindex;
							PsimagLite::Matrix<ComplexOrRealType> e(
							    nrow_e, ncol_e);

							den_kron_submatrix(nrow_a,
							                   ncol_a,
							                   a,
							                   nrow_b,
							                   ncol_b,
							                   b,
							                   nrindex,
							                   ncindex,
							                   rindex,
							                   cindex,
							                   e);

							/*
							 * --------------------------
							 * check E and D are the same
							 * --------------------------
							 */
							for (int je = 0; je < ncol_e; ++je) {
								for (int ie = 0; ie < nrow_e;
								     ++ie) {
									ComplexOrRealType eij
									    = e(ie, je);
									ComplexOrRealType dij
									    = d(ie, je);

									if (eij != dij) {
										INFO("nrow_A="
										     << nrow_a
										     << " ncol_A="
										     << ncol_a
										     << " nrow_B="
										     << nrow_b
										     << " ncol_B="
										     << ncol_b
										     << "  nrindex "
										     << nrindex
										     << " ncindex "
										     << ncindex
										     << '\n'
										     << " ie " << ie
										     << " je " << je
										     << " eij "
										     << eij
										     << " dij "
										     << dij);
										REQUIRE(eij == dij);
									}
								}
							}

							/*
							 * ------------------------------
							 * form sparse version of E from
							 * sparse version of A, B
							 * ------------------------------
							 */

							PsimagLite::CrsMatrix<ComplexOrRealType> e(
							    nrow_e, ncol_e);
							const int max_nnz_e = 1 + den_nnz(e);
							csr_kron_submatrix(a,
							                   b,
							                   nrindex,
							                   ncindex,
							                   max_nnz_e,
							                   rindex,
							                   cindex,
							                   e);

							/*
							 * ---------------------------------
							 * convert from sparse back to dense
							 * ---------------------------------
							 */
							PsimagLite::Matrix<ComplexOrRealType> se(
							    nrow_e, ncol_e);

							crsMatrixToFullMatrix(se, e);

							/*
							 * -----------------------------------------
							 * check E(ie,je) and SE(ie,je) are the same
							 * -----------------------------------------
							 */

							for (int je = 0; je < ncol_e; ++je) {
								for (int ie = 0; ie < nrow_e;
								     ++ie) {
									RealType diff = std::abs(
									    e(ie, je) - se(ie, je));
									const RealType tol = 1.0
									    / (1000.0 * 1000.0);
									if (diff > tol) {
										INFO("nrow_E="
										     << nrow_e
										     << " ncol_E="
										     << ncol_e

										     << " E(" << ie
										     << "," << je
										     << " "
										     << e(ie, je)
										     << " SE(" << ie
										     << "," << je
										     << ")"
										     << se(ie, je));
										REQUIRE(diff
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

TEST_CASE("kron_submatrix_test2_double", "[kron][submatrix]") { runKronSubmatrixChecks<double>(); }
TEST_CASE("kron_submatrix_test2_complex", "[kron][submatrix]")
{
	runKronSubmatrixChecks<std::complex<double>>();
}

int main(int argc, char* argv[])
{
	Kokkos::initialize(argc, argv);
	int result = Catch::Session().run(argc, argv);
	Kokkos::finalize();
	return result;
}
