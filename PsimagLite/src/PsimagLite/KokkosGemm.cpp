#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <complex>
#include <stdexcept>
#include <type_traits>

#include "KokkosGemm.h"
#include "KokkosType.h"

template <typename Scalar, typename IntegerForBlasType>
inline void PsimagLite::kokkosGemm(char               transa,
                                   char               transb,
                                   IntegerForBlasType m,
                                   IntegerForBlasType n,
                                   IntegerForBlasType k,
                                   const Scalar&      alpha,
                                   const Scalar*      A,
                                   IntegerForBlasType lda,
                                   const Scalar*      B,
                                   IntegerForBlasType ldb,
                                   const Scalar&      beta,
                                   Scalar*            C,
                                   IntegerForBlasType ldc)
{
	Kokkos::Profiling::ScopedRegion scoped_region("PsimagLite::kokkosGemm");
	int                             m       = static_cast<int>(m);
	int                             n       = static_cast<int>(n);
	int                             k       = static_cast<int>(k);
	int                             lda_val = static_cast<int>(lda);
	int                             ldb_val = static_cast<int>(ldb);
	int                             ldc_val = static_cast<int>(ldc);

	// Normalize trans flags
	char ta = transa ? transa : 'N';
	char tb = transb ? transb : 'N';
	if (ta >= 'a' && ta <= 'z')
		ta = char(ta - 'a' + 'A');
	if (tb >= 'a' && tb <= 'z')
		tb = char(tb - 'a' + 'A');

	int req_lda = (ta == 'N') ? std::max(1, m) : std::max(1, k);
	int req_ldb = (tb == 'N') ? std::max(1, k) : std::max(1, n);
	int req_ldc = std::max(1, m);
	if (lda_val < req_lda || ldb_val < req_ldb || ldc_val < req_ldc) {
		throw std::runtime_error("kokkosGemm: invalid leading dimension");
	}

	// Determine Kokkos scalar type
	using KokkosScalar = KokkosType<Scalar>::type;

	Kokkos::DefaultExecutionSpace exec;
	decltype(exec)::memory_space  mem;

	// allow padded leading dimensions (ldaVal/ldbVal/ldcVal >= required)
	if (lda_val < req_lda || ldb_val < req_ldb || ldc_val < req_ldc) {
		throw std::runtime_error("kokkosGemm: invalid leading dimension");
	}

	using Pair = Kokkos::pair<int, int>;

	// Create host unmanaged views that reflect the actual storage (use lda/ldb as the first
	// extent) and create subviews representing the logical matrix sizes.
	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	     aview_op(reinterpret_cast<const KokkosScalar*>(A), lda_val, (ta == 'N' ? k : m));
	auto aop_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), aview_op);
	auto aop_logical = (ta == 'N') ? Kokkos::subview(aop_device, Pair(0, m), Pair(0, k))
	                               : Kokkos::subview(aop_device, Pair(0, k), Pair(0, m));

	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	     bview_op(reinterpret_cast<const KokkosScalar*>(B), ldb_val, (tb == 'N' ? n : k));
	auto bop_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), bview_op);
	auto bop_logical = (tb == 'N') ? Kokkos::subview(bop_device, Pair(0, k), Pair(0, n))
	                               : Kokkos::subview(bop_device, Pair(0, n), Pair(0, k));

	// Create C view that reflects storage with possible padding (ldcVal >= M)
	Kokkos::View<KokkosScalar**, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
	     cview(reinterpret_cast<KokkosScalar*>(C), ldc_val, n);
	auto cop_device = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), cview);
	auto cop_logical = Kokkos::subview(cop_device, Pair(0, m), Pair(0, n));

	const char trans_a2[2] = { ta, '\0' };
	const char trans_b2[2] = { tb, '\0' };
	KokkosBlas::gemm(
	    exec, trans_a2, trans_b2, alpha, aop_logical, bop_logical, beta, cop_logical);
	Kokkos::deep_copy(exec, cview, cop_device);
	exec.fence();
}

#define PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(SCALAR, INTEGER)                                        \
	template void PsimagLite::kokkosGemm(char          transa,                                 \
	                                     char          transb,                                 \
	                                     INTEGER       m,                                      \
	                                     INTEGER       n,                                      \
	                                     INTEGER       k,                                      \
	                                     const SCALAR& alpha,                                  \
	                                     const SCALAR* A,                                      \
	                                     INTEGER       lda,                                    \
	                                     const SCALAR* B,                                      \
	                                     INTEGER       ldb,                                    \
	                                     const SCALAR& beta,                                   \
	                                     SCALAR*       C,                                      \
	                                     INTEGER       ldc)

#ifndef PSI_BLAS_64
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(double, int);
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(std::complex<double>, int);
#else
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(double, long int);
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(std::complex<double>, long int);
#endif

#undef PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM
