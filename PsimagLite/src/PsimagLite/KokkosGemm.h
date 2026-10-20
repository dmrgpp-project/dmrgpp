#ifndef PSIMAG_KOKKOS_GEMM_H
#define PSIMAG_KOKKOS_GEMM_H

#include <Kokkos_Core.hpp>

namespace PsimagLite {

template <typename Scalar, typename KokkosScalar, typename IntegerForBlasType>
void kokkos_gemm(const Kokkos::DefaultExecutionSpace& exec,
                                    char                                 transa,
                                    char                                 transb,
                                    IntegerForBlasType                   m,
                                    IntegerForBlasType                   n,
                                    IntegerForBlasType                   k,
                                    const Scalar&                        alpha,
                                    Kokkos::View<const KokkosScalar*, Kokkos::DefaultExecutionSpace> A,
                                    IntegerForBlasType                   lda,
                                    Kokkos::View<const KokkosScalar*, Kokkos::DefaultExecutionSpace> B,
                                    IntegerForBlasType                   ldb,
                                    const Scalar&                        beta,
                                    Kokkos::View<KokkosScalar*, Kokkos::DefaultExecutionSpace> C,
                                    IntegerForBlasType                   ldc);

template <typename Scalar, typename IntegerForBlasType>
void kokkos_gemm(const Kokkos::DefaultExecutionSpace& exec,
                 char                                 transa,
                 char                                 transb,
                 IntegerForBlasType                   m,
                 IntegerForBlasType                   n,
                 IntegerForBlasType                   k,
                 const Scalar&                        alpha,
                 const Scalar*                        A,
                 IntegerForBlasType                   lda,
                 const Scalar*                        B,
                 IntegerForBlasType                   ldb,
                 const Scalar&                        beta,
                 Scalar*                              C,
                 IntegerForBlasType                   ldc);
} // namespace PsimagLite

#endif // PSIMAG_KOKKOS_GEMM_H
