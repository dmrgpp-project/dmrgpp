#ifndef BATCHED_GEMM_KOKKOS_H
#define BATCHED_GEMM_KOKKOS_H
// Don't include this file directly; use BatchedGemmInclude.hh

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/KokkosType.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_TeamVector_Impl.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <vector>

namespace Dmrg {

/*
 * Batched GEMM using KokkosKernels with a sparse-batch layout.
 *
 * Algorithm (mirrors apply_Htarget_sparse_not_ready.cpp):
 *
 *  Pass 1: For each non-zero triple (ipatch, jpatch, ioperator):
 *     BXbatch[ip][:, colBX:colBX+lps[jp]] =
 *         Bbatch[ip][:, colB:colB+rps[jp]] * X[jp]
 *
 *  Pass 2: For each ipatch:
 *     Y[ip] = BXbatch[ip] * Abatch[ip]^T
 *
 * Abatch[ip] = left-operator (xc) matrices packed column-by-column (ldA padded).
 * Bbatch[ip] = right-operator (yc) matrices packed column-by-column (ldB padded).
 * BXbatch[ip] is the intermediate work buffer (ldB padded, same col-width as Abatch).
 *
 * All Abatch and Bbatch data live on the device for the lifetime of the object.
 * BXbatch, vin, vout are (re-)written on each matrixVector call.
 */

template <typename InitKronType> class BatchedGemmKokkos {

	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using GenIjPatchType          = typename InitKronType::GenIjPatchType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using VectorType              = typename MatrixDenseOrSparseType::VectorType;
	using ComplexOrRealType       = typename VectorType::value_type;
	using MatrixType              = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorMatrixType = typename PsimagLite::Vector<std::unique_ptr<MatrixType>>::Type;
	using VectorSizeType   = PsimagLite::Vector<SizeType>::Type;

	// Kokkos types
	using KokkosScalar   = typename PsimagLite::KokkosType<ComplexOrRealType>::type;
	using ExecutionSpace = Kokkos::DefaultExecutionSpace;
	using MemorySpace    = typename ExecutionSpace::memory_space;
	using ScalarView     = Kokkos::View<KokkosScalar*, MemorySpace>;

	// Needs to be public so we can use lambdas with Cuda.
public:

	// Compact struct holding all parameters for one batched GEMM call.
	// Stored in device memory and accessed inside the kernel.
	struct GemmArgs {
		int       m, n, k; // GEMM dimensions
		int       lda, ldb, ldc; // leading dimensions (column-major strides)
		long long a_off, b_off, c_off; // element offsets into device flat arrays
	};

private:

	using DevArgsView = Kokkos::View<GemmArgs*, MemorySpace>;

	inline static constexpr int ialign_ = 32;

	// Pass 2 GEMMs have very few batches (one per patch, == npatches, which can be as low as
	// single digits) but each has a very long contraction (k) dimension, so a single GEMM per
	// patch leaves the GPU with far too few concurrent thread-blocks to reach good occupancy.
	// We split the k dimension of each patch's GEMM into independent slices and accumulate
	// their partial results directly into d_vout_ via atomic adds (see matrixVector()). This
	// multiplies the number of pass-2 batches (and thus concurrent thread-blocks) without
	// changing the per-GEMM output shape (so the favorable Blocked-algo memory access pattern
	// is preserved).
	//
	// Splitting into a *fixed count* of chunks per patch under-performs when patches have very
	// different k (contraction) sizes: every patch still gets the same number of chunks, so
	// patches with large k get large (slow) chunks while patches with small k get tiny (fast)
	// chunks. Since all chunks across all patches are launched in the same kernel, the
	// wall-clock time is dominated by the few largest chunks. Instead we split by a *fixed
	// column width* (kPass2ColChunk_): each patch gets ceil(k_ip / kPass2ColChunk_) chunks, so
	// patches with larger k automatically get proportionally more (similarly-sized) chunks.
	// This balances per-team work across the whole grid regardless of how unevenly k is
	// distributed across patches.
	//
	// Each chunk-team computes its partial result into a small *per-team* scratch buffer backed
	// by Kokkos' global- memory scratch pool (level 1), whose total size is bounded by the
	// number of *concurrently resident* teams (a small, GPU-occupancy-bounded constant) times
	// the per-team scratch size -- NOT by the total number of chunks across the whole batch.
	// The team then atomically adds its scratch buffer into d_vout_ (zeroed before Pass 2).
	//
	// However, a single patch's own output (m x n = rps[ip] * lps[ip]) can itself be very large
	// for big systems (observed: a single patch's scratch requirement reached tens of GiB when
	// multiplied across the GPU's many concurrently-resident teams). To keep the per-team
	// scratch size bounded by a small constant regardless of problem size, we additionally tile
	// each patch's (m x n) output into at most kPass2TileDim_ x kPass2TileDim_ blocks (in
	// addition to the k-chunking above); output tiles/rows are independent (no reduction needed
	// across tiles, only across k-chunks of the same tile), so this is a straightforward extra
	// split, analogous to Pass 1's row-splitting.

	// The exact values have been tuned for a Grace Hopper 200 and the 345 test case. The
	// bench_BatchedGemm executable can be used for tuning the values for different
	// architectures if necessary. It turned out that they were also suitable for a MI300A.
	inline static constexpr int kPass2ColChunk_ = 32;
	inline static constexpr int kPass2TileDim_  = 64;

	// Pass 1 batches one GEMM per non-zero (ip, jp, k) connection triple, but the number of
	// such triples (tens to low thousands, see setup_) is often far smaller than what's needed
	// to fill a modern GPU with concurrent thread-blocks, even though many of those GEMMs have
	// a sizeable number of output rows (m = rps[ip], which can range into the hundreds). Since
	// output rows are fully independent (no reduction needed, unlike the k-split used for Pass
	// 2), we split each triple's m dimension into row-blocks of at most kPass1RowChunk_ rows
	// and emit one GEMM per row-block, multiplying the number of pass-1 teams for large-m
	// triples while leaving small-m triples (m <= kPass1RowChunk_) unsplit.
	inline static constexpr int kPass1RowChunk_ = 16;

	// By default (and when using Kokkos::AUTO), Kokkos uses 1 as vector_length which means that
	// all the GEMMs would be executed by only one hardware thread. Just as for the other
	// parameters above, it turned out that we get get better performance by letting multiple
	// hardware threads collaborate on a GEMM.
	inline static constexpr int kVectorLength_ = 8;

	// The work is partitioned across all available devices: patches are split
	// into contiguous ranges, balanced by an approximate per-patch cost, and
	// each range gets its own DeviceShard below.
	//
	// The input vector vin is comparatively small, so it is simply
	// replicated in full to every shard (Pass 1 of any patch ip may
	// reference any other patch jp's slice of vin, so a shard cannot know in
	// advance which slices it will need). Abatch/Bbatch/BXbatch/vout, which
	// can be much larger, are only ever built/sized for a shard's own
	// patches. Because patches are assigned to exactly one shard each, every
	// shard's vout contribution lands in a disjoint region of the final
	// output vector, so no cross-device reduction step is needed.
	struct DeviceShard {
		ExecutionSpace exec;
		SizeType       patchBegin = 0, patchEnd = 0; // half-open patch range
		SizeType voutHostOffset = 0; // offset of this shard's slice in the host vout array
		SizeType voutLocalSize  = 0; // size of this shard's slice
		SizeType nbatch1        = 0; // number of pass-1 GEMMs in this shard
		SizeType nbatch2        = 0; // number of pass-2 GEMMs in this shard

		ScalarView         d_flatAbatch; // left-operator matrices, device, persistent
		ScalarView         d_flatBbatch; // right-operator matrices, device, persistent
		mutable ScalarView d_flatBXbatch; // intermediate BX work buffer, per-call
		mutable ScalarView d_vin; // FULL input vector replicated on this device, per-call
		mutable ScalarView d_vout; // this shard's local-indexed output slice, per-call
		DevArgsView        d_pass1; // pass-1 GEMM parameters, persistent
		DevArgsView        d_pass2; // pass-2 GEMM parameters, persistent

		// SYCL doesn't need to manually manage the lifecycle of sycl::queues
#if defined(KOKKOS_ENABLE_CUDA)
		cudaStream_t stream = nullptr;
#endif
#if defined(KOKKOS_ENABLE_HIP)
		hipStream_t stream = nullptr;
#endif
	};

public:

	BatchedGemmKokkos(const InitKronType& initKron)
	    : initKron_(initKron)
	    , progress_("BatchedGemmKokkos")
	{
		if (!enabled())
			return;
		setup_();
	}

	~BatchedGemmKokkos()
	{
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
		for (auto& shard : shards_) {
			if (shard.stream != nullptr) {
				shard.exec.fence("BatchedGemmKokkos::~BatchedGemmKokkos");
#if defined(KOKKOS_ENABLE_HIP)
				(void)hipStreamDestroy(shard.stream);
#else
				(void)cudaStreamDestroy(shard.stream);
#endif
				shard.stream = nullptr;
			}
		}
#endif
	}

	bool enabled() const { return initKron_.params().options.isSet("BatchedGemm"); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		assert(enabled());

		Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::matrixVector");

		const SizeType totalXY = vin.size();
		assert(vout.size() == totalXY);
		assert(!shards_.empty());
		assert(shards_[0].d_vin.extent(0) == totalXY);

		using HVin
		    = Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		HVin hvin(reinterpret_cast<const KokkosScalar*>(vin.data()), totalXY);

		// --- Launch phase -----------------------------------------------------
		// Enqueue the H2D copy, Pass 1 and Pass 2 for every shard WITHOUT
		// fencing in between: each shard uses its own execution-space
		// instance (bound to its own device/stream), so leaving all of them
		// unfenced here lets independent shards' GPU work actually overlap
		// across physical devices. Fencing happens only in the completion
		// phase below, after every shard's work has been enqueued.
		for (const DeviceShard& shard : shards_) {
			const ExecutionSpace& exec = shard.exec;

			// --- H2D: replicate the FULL vin onto this shard's device ------
			Kokkos::deep_copy(exec, shard.d_vin, hvin);

			// Note: d_flatBXbatch is NOT explicitly zeroed here. Every
			// element is written with beta=0 by TeamVectorGemmInternal
			// inside Pass 1 (which zero-fills C whenever beta==0, even for
			// k==0 slices), so an upfront deep_copy zeroing would be
			// redundant extra global-memory traffic. d_vout, however, IS
			// explicitly zeroed just before Pass 2 below, since Pass 2's
			// split-k chunks accumulate into it via atomic adds (see
			// kPass2ColChunk_ comment above) rather than each writing a
			// disjoint region.

			// --- Pass 1: BXbatch[ip] = Bbatch[ip] * X[jp]  (NoT x NoT) -----
			{
				const ScalarView  flatBbatch  = shard.d_flatBbatch;
				const ScalarView  vin_dev     = shard.d_vin;
				const ScalarView  flatBXbatch = shard.d_flatBXbatch;
				const DevArgsView args        = shard.d_pass1;

				using MemberType =
				    typename Kokkos::TeamPolicy<ExecutionSpace>::member_type;
				Kokkos::parallel_for(
				    "BatchedGemmKokkos_Pass1",
				    Kokkos::TeamPolicy<ExecutionSpace>(
				        exec,
				        static_cast<int>(shard.nbatch1),
				        Kokkos::AUTO,
				        kVectorLength_),
				    KOKKOS_LAMBDA(const MemberType& member) {
					    const int       i  = member.league_rank();
					    const GemmArgs& ag = args(i);

					    using UV = Kokkos::View<KokkosScalar**,
					                            Kokkos::LayoutLeft,
					                            MemorySpace,
					                            Kokkos::MemoryUnmanaged>;

					    // A = Bbatch[ip] block: (lda, k) padded -> subview (m,
					    // k) active rows
					    UV   A_full(flatBbatch.data() + ag.a_off, ag.lda, ag.k);
					    auto A = Kokkos::subview(
					        A_full, Kokkos::make_pair(0, ag.m), Kokkos::ALL());

					    // B = X[jp]: (k, n) exact -- no padding for X slice
					    UV B(vin_dev.data() + ag.b_off, ag.ldb, ag.n);

					    // C = BXbatch[ip] block: (ldc, n) padded -> subview (m,
					    // n) active rows
					    UV C_full(flatBXbatch.data() + ag.c_off, ag.ldc, ag.n);
					    auto C = Kokkos::subview(
					        C_full, Kokkos::make_pair(0, ag.m), Kokkos::ALL());

					    KokkosBatched::TeamVectorGemm<
					        MemberType,
					        KokkosBatched::Trans::NoTranspose,
					        KokkosBatched::Trans::NoTranspose,
					        KokkosBatched::Algo::Gemm::Blocked>::
					        invoke(member,
					               KokkosScalar(1),
					               A,
					               B,
					               KokkosScalar(0),
					               C);
				    });
			}

			// --- Pass 2: Y[ip] += BXbatch[ip]_chunk * Abatch[ip]_chunk^T
			// -------------- (split-k: each patch's contraction dim is divided into
			// independent column-width-bounded slices for better GPU occupancy;
			// see kPass2ColChunk_ comment above. Each chunk-team computes its
			// partial (m x n) result into a small per-team scratch buffer -- NOT
			// a global buffer sized by total chunk count -- then atomically
			// adds it into d_vout_, so d_vout_ must be zeroed first.)
			{
				Kokkos::deep_copy(exec, shard.d_vout, KokkosScalar(0));

				const ScalarView  flatBXbatch = shard.d_flatBXbatch;
				const ScalarView  flatAbatch  = shard.d_flatAbatch;
				const ScalarView  vout_dev    = shard.d_vout;
				const DevArgsView args        = shard.d_pass2;
				// Fixed, compile-time-bounded per-team scratch size: every
				// Pass 2 GEMM's output tile is at most kPass2TileDim_ x
				// kPass2TileDim_ (see setup_()), so this does NOT grow with
				// problem size.
				constexpr size_t scratchBytesPerTeam
				    = static_cast<size_t>(kPass2TileDim_) * kPass2TileDim_
				    * sizeof(KokkosScalar);

				using MemberType =
				    typename Kokkos::TeamPolicy<ExecutionSpace>::member_type;
				using ScratchView
				    = Kokkos::View<KokkosScalar**,
				                   Kokkos::LayoutLeft,
				                   typename ExecutionSpace::scratch_memory_space,
				                   Kokkos::MemoryUnmanaged>;
				Kokkos::parallel_for(
				    "BatchedGemmKokkos_Pass2",
				    Kokkos::TeamPolicy<ExecutionSpace>(
				        exec,
				        static_cast<int>(shard.nbatch2),
				        Kokkos::AUTO,
				        kVectorLength_)
				        .set_scratch_size(
				            1,
				            Kokkos::PerTeam(static_cast<int>(scratchBytesPerTeam))),
				    KOKKOS_LAMBDA(const MemberType& member) {
					    const int       i  = member.league_rank();
					    const GemmArgs& ag = args(i);

					    // Note: we deliberately do NOT special-case ag.k == 0
					    // here. TeamVectorGemmInternal zero-fills C whenever
					    // beta == 0, even if k == 0, so every element of this
					    // chunk's output scratch is always written by the
					    // invoke() call below.

					    using UV = Kokkos::View<KokkosScalar**,
					                            Kokkos::LayoutLeft,
					                            MemorySpace,
					                            Kokkos::MemoryUnmanaged>;

					    // A = BXbatch[ip] chunk: (lda, k) padded -> subview (m,
					    // k)
					    UV A_full(flatBXbatch.data() + ag.a_off, ag.lda, ag.k);
					    auto A = Kokkos::subview(
					        A_full, Kokkos::make_pair(0, ag.m), Kokkos::ALL());

					    // B = Abatch[ip] chunk: (ldb, k) padded -> subview (n,
					    // k); will be transposed
					    UV   B_full(flatAbatch.data() + ag.b_off, ag.ldb, ag.k);
					    auto B = Kokkos::subview(
					        B_full, Kokkos::make_pair(0, ag.n), Kokkos::ALL());

					    // C = this chunk's partial Y[ip]: small per-team
					    // scratch buffer (m x n, contiguous), NOT a slice of a
					    // global buffer sized by total chunk count -- see
					    // comment above.
					    ScratchView C(member.team_scratch(1), ag.m, ag.n);

					    // Pass 2 GEMMs have a very "thin" output (m, n are the
					    // per-patch sizes) with a very long contraction
					    // dimension k. The Blocked algorithm reuses loaded A/B
					    // tiles across the register-blocked inner loop, which
					    // is far more memory-bandwidth efficient for this
					    // thin/long GEMM shape than assigning one thread per
					    // output element (tried Unblocked: ~3x slower here
					    // because every thread re-reads all of A/B from global
					    // memory).
					    KokkosBatched::TeamVectorGemm<
					        MemberType,
					        KokkosBatched::Trans::NoTranspose,
					        KokkosBatched::Trans::Transpose,
					        KokkosBatched::Algo::Gemm::Blocked>::
					        invoke(member,
					               KokkosScalar(1),
					               A,
					               B,
					               KokkosScalar(0),
					               C);

					    // Accumulate this chunk's partial result into the
					    // shared global output via atomic adds (multiple chunks
					    // of the same patch race on the same d_vout_ region).
					    // Note: ag.ldc is the FULL patch's row count (rps[ip]),
					    // used as the addressing stride into d_vout_ -- NOT the
					    // same as ag.m, which is only this tile's row count
					    // (see kPass2TileDim_ comment above).
					    member.team_barrier();
					    Kokkos::parallel_for(
					        Kokkos::TeamThreadRange(member, ag.n),
					        [&](const int j)
					        {
						        Kokkos::parallel_for(
						            Kokkos::ThreadVectorRange(member, ag.m),
						            [&](const int r)
						            {
							            Kokkos::atomic_add(
							                &vout_dev(
							                    ag.c_off + r
							                    + static_cast<
							                          long long>(j)
							                        * ag.ldc),
							                C(r, j));
						            });
					        });
				    });
			}
		}

		// --- Completion phase ---------------------------------------------
		// Every shard's H2D copy, Pass 1 and Pass 2 have now been enqueued
		// (but not necessarily finished) above. Enqueue each shard's D2H
		// copy of its own disjoint vout slice next (again without fencing in
		// between, for the same cross-device overlap reason as above), and
		// only fence -- waiting for genuine completion -- once every
		// shard's work, including the D2H copy, has been submitted.
		for (const DeviceShard& shard : shards_) {
			using HVout = Kokkos::
			    View<KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			HVout hv(reinterpret_cast<KokkosScalar*>(vout.data())
			             + shard.voutHostOffset,
			         shard.voutLocalSize);
			Kokkos::deep_copy(shard.exec, hv, shard.d_vout);
		}
		for (const DeviceShard& shard : shards_)
			shard.exec.fence("BatchedGemmKokkos::matrixVector");
	}

private:

	static int iceil(int x, int n) { return (x + n - 1) / n; }

	// Dense matrix reference: pointer + source dimensions (no padding).
	// For sparse matrices, the expansion is heap-allocated and owned by garbage_.
	struct MatRef {
		ComplexOrRealType* ptr;
		int                rows; // = source leading dimension
		int                cols;
	};

	MatRef getMatRef(const MatrixDenseOrSparseType& mat)
	{
		if (mat.isZero())
			return { nullptr, 0, 0 };
		if (!mat.isDense()) {
			MatrixType* dense = new MatrixType();
			crsMatrixToFullMatrix(*dense, mat.sparse());
			garbage_.emplace_back(dense);
			return { &(*dense)(0, 0),
				 static_cast<int>(dense->rows()),
				 static_cast<int>(dense->cols()) };
		}
		const MatrixType& d = mat.dense();
		return { const_cast<ComplexOrRealType*>(&d(0, 0)),
			 static_cast<int>(d.rows()),
			 static_cast<int>(d.cols()) };
	}

	// Needs to be public so we can use lambdas with Cuda.
public:

	struct PackTask {
		const ComplexOrRealType* src;
		ComplexOrRealType*       dst;
		int                      m, n, lds, ldd;
	};

	void setup_()
	{
		Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup");

		const SizeType npatches  = initKron_.numberOfPatches(InitKronType::OLD);
		const SizeType noperator = initKron_.connections();

		// Per-patch sizes and padded leading dimensions.
		VectorSizeType lps(npatches, 0); // left patch size
		VectorSizeType rps(npatches, 0); // right patch size
		VectorSizeType ldA(npatches, 0); // padded lps -> Abatch leading dim
		VectorSizeType ldB(npatches, 0); // padded rps -> Bbatch / BXbatch leading dim
		VectorSizeType xyStart(npatches + 1, 0); // patch start offset in vin/vout

		{
			Kokkos::Profiling::ScopedRegion region(
			    "BatchedGemmKokkos::setup::patch_sizes");

			for (SizeType ip = 0; ip < npatches; ++ip) {
				const SizeType lg
				    = initKron_.patch(InitKronType::NEW, GenIjPatchType::LEFT)[ip];
				const SizeType rg
				    = initKron_.patch(InitKronType::NEW, GenIjPatchType::RIGHT)[ip];
				const int L1
				    = initKron_.lrs(InitKronType::NEW).left().partition(lg);
				const int L2
				    = initKron_.lrs(InitKronType::NEW).left().partition(lg + 1);
				const int R1
				    = initKron_.lrs(InitKronType::NEW).right().partition(rg);
				const int R2
				    = initKron_.lrs(InitKronType::NEW).right().partition(rg + 1);
				lps[ip] = static_cast<SizeType>(L2 - L1);
				rps[ip] = static_cast<SizeType>(R2 - R1);
				ldA[ip] = static_cast<SizeType>(
				    ialign_ * iceil(static_cast<int>(lps[ip]), ialign_));
				ldB[ip] = static_cast<SizeType>(
				    ialign_ * iceil(static_cast<int>(rps[ip]), ialign_));
			}

			xyStart[0] = 0;
			for (SizeType ip = 0; ip < npatches; ++ip)
				xyStart[ip + 1] = xyStart[ip] + lps[ip] * rps[ip];
		}
		// Column widths:
		//   AbatchCols[ip] = sum of lps[jp] over all non-zero (jp, k) connections for ip
		//   BbatchCols[ip] = sum of rps[jp] over all non-zero (jp, k) connections for ip
		VectorSizeType AbatchCols(npatches, 0);
		VectorSizeType BbatchCols(npatches, 0);

		{
			Kokkos::Profiling::ScopedRegion region(
			    "BatchedGemmKokkos::setup::patch_pointers");
			for (SizeType ip = 0; ip < npatches; ++ip) {
				for (SizeType jp = 0; jp < npatches; ++jp) {
					for (SizeType k = 0; k < noperator; ++k) {
						const MatrixDenseOrSparseType* Ap
						    = initKron_.xc(k)(ip, jp);
						const MatrixDenseOrSparseType* Bp
						    = initKron_.yc(k)(ip, jp);
						if (!Ap || !Bp)
							continue;
						if (Ap->isZero() || Bp->isZero())
							continue;
						AbatchCols[ip] += lps[jp];
						BbatchCols[ip] += rps[jp];
					}
				}
			}
		}
		Kokkos::Profiling::pushRegion("BatchedGemmKokkos::setup::intermediate");
		// Flat offsets: each patch ip occupies a contiguous slice.
		VectorSizeType AbatchOff(npatches + 1, 0);
		VectorSizeType BbatchOff(npatches + 1, 0);
		VectorSizeType BXbatchOff(npatches + 1, 0); // ldB rows x AbatchCols columns

		for (SizeType ip = 0; ip < npatches; ++ip) {
			AbatchOff[ip + 1]  = AbatchOff[ip] + ldA[ip] * AbatchCols[ip];
			BbatchOff[ip + 1]  = BbatchOff[ip] + ldB[ip] * BbatchCols[ip];
			BXbatchOff[ip + 1] = BXbatchOff[ip] + ldB[ip] * AbatchCols[ip];
		}

		const SizeType totalAbatch  = AbatchOff[npatches];
		const SizeType totalBbatch  = BbatchOff[npatches];
		const SizeType totalBXbatch = BXbatchOff[npatches];

		// Pass-2 split-k chunk count: patch ip's contraction dim
		// (AbatchCols[ip]) is split into ceil(AbatchCols[ip] /
		// kPass2ColChunk_) chunks (at least 1), so patches with a larger
		// contraction dimension automatically get proportionally more
		// (similarly-sized) chunks -- see kPass2ColChunk_ comment above.
		// All chunks of the same patch/tile atomically accumulate into the
		// same (m x n) region of d_vout_ (see matrixVector()).
		VectorSizeType voutChunkCount(npatches, 0);
		for (SizeType ip = 0; ip < npatches; ++ip) {
			voutChunkCount[ip] = static_cast<SizeType>(
			    std::max(1, iceil(static_cast<int>(AbatchCols[ip]), kPass2ColChunk_)));
		}

		Kokkos::Profiling::popRegion();

		// --- Partition patches across available GPU devices --------------
		// Contiguous patch ranges, balanced by an approximate per-patch cost
		// (AbatchCols[ip]*lps[ip]*rps[ip], a proxy for Pass 2's flop count --
		// the phase profiling found to dominate GPU time), so that each
		// shard below gets roughly equal work. See the DeviceShard comment
		// for why a contiguous-range partition (rather than e.g. a
		// round-robin one) lets each shard's vout slice be copied back with
		// a single contiguous D2H copy and avoids any cross-device
		// reduction.
		std::vector<std::pair<SizeType, SizeType>> shardRanges;
		{
			const int numDevices
			    = std::max(1,
			               std::min(Kokkos::num_devices(),
			                        static_cast<int>(std::max<SizeType>(1, npatches))));
			VectorSizeType      bound(numDevices + 1, 0);
			std::vector<double> cost(std::max<SizeType>(1, npatches), 0.0);
			double              total = 0.0;
			for (SizeType ip = 0; ip < npatches; ++ip) {
				cost[ip]
				    = static_cast<double>(AbatchCols[ip]) * lps[ip] * rps[ip] + 1.0;
				total += cost[ip];
			}
			SizeType ip  = 0;
			double   acc = 0.0;
			for (int s = 0; s < numDevices; ++s) {
				if (s + 1 == numDevices) {
					bound[s + 1] = npatches;
					continue;
				}
				const double target = total * (s + 1) / numDevices;
				while (ip < npatches && acc < target) {
					acc += cost[ip];
					++ip;
				}
				bound[s + 1] = ip;
			}
			for (int s = 0; s < numDevices; ++s) {
				if (bound[s + 1] > bound[s])
					shardRanges.emplace_back(bound[s], bound[s + 1]);
			}
			if (shardRanges.empty())
				shardRanges.emplace_back(SizeType(0), npatches);
		}

		const SizeType numShards = shardRanges.size();
		shards_.clear();
		shards_.resize(numShards);
		for (SizeType s = 0; s < numShards; ++s) {
			shards_[s].patchBegin = shardRanges[s].first;
			shards_[s].patchEnd   = shardRanges[s].second;
#if defined(KOKKOS_ENABLE_HIP)
			if (numShards > 1 && hipSetDevice(static_cast<int>(s)) == hipSuccess) {
				hipStream_t stream = nullptr;
				if (hipStreamCreate(&stream) == hipSuccess) {
					shards_[s].stream = stream;
					shards_[s].exec   = ExecutionSpace(stream);
				}
			}
#endif
#if defined(KOKKOS_ENABLE_CUDA)
			if (numShards > 1 && cudaSetDevice(static_cast<int>(s)) == cudaSuccess) {
				cudaStream_t stream = nullptr;
				if (cudaStreamCreate(&stream) == cudaSuccess) {
					shards_[s].stream = stream;
					shards_[s].exec   = ExecutionSpace(stream);
				}
			}
#endif
		}

		using HostMirrorView = typename ScalarView::host_mirror_type;
		std::vector<HostMirrorView> hostAbatch(numShards);
		std::vector<HostMirrorView> hostBbatch(numShards);
		std::vector<PackTask>       packTasks; // combined across all shards
		packTasks.reserve(npatches * npatches);

		{
			Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup::gemm");

			for (SizeType s = 0; s < numShards; ++s) {
				DeviceShard&   shard      = shards_[s];
				const SizeType patchBegin = shard.patchBegin;
				const SizeType patchEnd   = shard.patchEnd;
				const SizeType nLocal     = patchEnd - patchBegin;

				// Local (shard-relative) flat offsets: patch ip (global)
				// maps to local index ip - patchBegin. Abatch/Bbatch/BXbatch
				// are only ever built for this shard's own patches.
				VectorSizeType AbatchOff(nLocal + 1, 0);
				VectorSizeType BbatchOff(nLocal + 1, 0);
				VectorSizeType BXbatchOff(nLocal + 1, 0);
				for (SizeType li = 0; li < nLocal; ++li) {
					const SizeType ip = patchBegin + li;
					AbatchOff[li + 1]
					    = AbatchOff[li] + ldA[ip] * AbatchCols[ip];
					BbatchOff[li + 1]
					    = BbatchOff[li] + ldB[ip] * BbatchCols[ip];
					BXbatchOff[li + 1]
					    = BXbatchOff[li] + ldB[ip] * AbatchCols[ip];
				}
				const SizeType totalAbatchShard  = AbatchOff[nLocal];
				const SizeType totalBbatchShard  = BbatchOff[nLocal];
				const SizeType totalBXbatchShard = BXbatchOff[nLocal];

				shard.d_flatAbatch = ScalarView(
				    Kokkos::view_alloc(
				        shard.exec, Kokkos::WithoutInitializing, "d_flatAbatch"),
				    totalAbatchShard);
				shard.d_flatBbatch = ScalarView(
				    Kokkos::view_alloc(
				        shard.exec, Kokkos::WithoutInitializing, "d_flatBbatch"),
				    totalBbatchShard);
				shard.d_flatBXbatch = ScalarView(
				    Kokkos::view_alloc(
				        shard.exec, Kokkos::WithoutInitializing, "d_flatBXbatch"),
				    totalBXbatchShard);

				hostAbatch[s] = Kokkos::create_mirror_view(
				    Kokkos::view_alloc(Kokkos::WithoutInitializing),
				    shard.d_flatAbatch);
				hostBbatch[s] = Kokkos::create_mirror_view(
				    Kokkos::view_alloc(Kokkos::WithoutInitializing),
				    shard.d_flatBbatch);
				auto& h_flatAbatch = hostAbatch[s];
				auto& h_flatBbatch = hostBbatch[s];

				std::vector<GemmArgs> pass1_args;
				pass1_args.reserve(nLocal * npatches * noperator);
				std::vector<GemmArgs> pass2_args;
				pass2_args.reserve(nLocal);

				for (SizeType ip = patchBegin; ip < patchEnd; ++ip) {
					const SizeType li = ip - patchBegin; // local patch index
					long long      colA
					    = 0; // column cursor in Abatch[ip] (= BXbatch[ip])
					long long colB = 0; // column cursor in Bbatch[ip]

					for (SizeType jp = 0; jp < npatches; ++jp) {
						for (SizeType k = 0; k < noperator; ++k) {
							const MatrixDenseOrSparseType* Amat
							    = initKron_.xc(k)(ip, jp);
							const MatrixDenseOrSparseType* Bmat
							    = initKron_.yc(k)(ip, jp);
							if (!Amat || !Bmat)
								continue;
							if (Amat->isZero() || Bmat->isZero())
								continue;

							MatRef Aref = getMatRef(
							    *Amat); // left operator (xc)
							MatRef Bref = getMatRef(
							    *Bmat); // right operator (yc)
							assert(Aref.ptr && Bref.ptr);

							const int mA  = static_cast<int>(lps[ip]);
							const int nAk = static_cast<int>(
							    lps[jp]); // A cols = BX cols per op
							const int mB  = static_cast<int>(rps[ip]);
							const int nBk = static_cast<int>(
							    rps[jp]); // B cols = X rows

							// Record pack task for A (left operator)
							PackTask tA;
							tA.src = Aref.ptr;
							tA.dst
							    = reinterpret_cast<ComplexOrRealType*>(
							          h_flatAbatch.data())
							    + AbatchOff[li]
							    + colA
							        * static_cast<long long>(ldA[ip]);
							tA.m   = mA;
							tA.n   = nAk;
							tA.lds = Aref.rows;
							tA.ldd = static_cast<int>(ldA[ip]);
							packTasks.push_back(tA);

							// Record pack task for B (right operator)
							PackTask tB;
							tB.src = Bref.ptr;
							tB.dst
							    = reinterpret_cast<ComplexOrRealType*>(
							          h_flatBbatch.data())
							    + BbatchOff[li]
							    + colB
							        * static_cast<long long>(ldB[ip]);
							tB.m   = mB;
							tB.n   = nBk;
							tB.lds = Bref.rows;
							tB.ldd = static_cast<int>(ldB[ip]);
							packTasks.push_back(tB);

							// Pass 1 GEMM: C(mB, nAk) =
							// Bbatch_block(mB, nBk) * X_jp(nBk, nAk)
							GemmArgs a1;
							a1.m   = mB;
							a1.n   = nAk;
							a1.k   = nBk;
							a1.lda = static_cast<int>(
							    ldB[ip]); // Bbatch leading dim
							a1.ldb
							    = nBk; // X[jp] no padding: ld = rps[jp]
							a1.ldc = static_cast<int>(
							    ldB[ip]); // BXbatch leading dim
							a1.a_off
							    = static_cast<long long>(BbatchOff[li])
							    + colB
							        * static_cast<long long>(ldB[ip]);
							// b_off indexes into d_vin, which is a FULL
							// replica of vin on every shard's device,
							// so it uses the GLOBAL xyStart (not
							// shard-local).
							a1.b_off
							    = static_cast<long long>(xyStart[jp]);
							a1.c_off
							    = static_cast<long long>(BXbatchOff[li])
							    + colA
							        * static_cast<long long>(ldB[ip]);

							// Split the m (=mB) dimension into
							// independent row-blocks: rows of C/A are
							// not a contraction dimension, so each
							// row-block is a fully self-contained GEMM
							// and no reduction step is needed (unlike
							// the k-split used for Pass 2).
							for (int r0 = 0; r0 < a1.m;
							     r0 += kPass1RowChunk_) {
								const int mc = std::min(
								    kPass1RowChunk_, a1.m - r0);
								GemmArgs ar = a1;
								ar.m        = mc;
								ar.a_off    = a1.a_off + r0;
								ar.c_off    = a1.c_off + r0;
								pass1_args.push_back(ar);
							}

							colA += nAk;
							colB += nBk;
						}
					}

					// Pass 2 GEMM: Y[ip](mB, mA) = BXbatch[ip](mB, k) *
					// Abatch[ip](mA, k)^T
					const int totalCols
					    = static_cast<int>(colA); // = AbatchCols[ip]
					GemmArgs a2;
					a2.m   = static_cast<int>(rps[ip]);
					a2.n   = static_cast<int>(lps[ip]);
					a2.k   = totalCols;
					a2.lda = static_cast<int>(ldB[ip]); // BXbatch leading dim
					a2.ldb = static_cast<int>(ldA[ip]); // Abatch leading dim
					a2.ldc = static_cast<int>(
					    rps[ip]); // Y[ip] no padding: ld = rps[ip]

					// Split the k (=totalCols) dimension into
					// voutChunkCount[ip] slices (proportional to this patch's
					// own k -- see kPass2ColChunk_ comment above), AND tile the
					// (m x n) output into at most kPass2TileDim_ x
					// kPass2TileDim_ blocks (see kPass2TileDim_ comment above,
					// bounds the Pass 2 per-team scratch buffer size). All
					// k-slices of the same (patch, output tile) write to the
					// SAME (mTile x nTile) region of this shard's d_vout via
					// atomic add in matrixVector(); output tiles are
					// independent (no reduction needed across tiles). Slices
					// with 0 k-width still get emitted so that
					// TeamVectorGemmInternal's beta==0 zero-fill produces a
					// well-defined (zero) partial result to add.
					const SizeType nChunks = voutChunkCount[ip];
					const int      chunkCols
					    = iceil(totalCols,
					            static_cast<int>(nChunks)); // columns per chunk
					for (int r0 = 0; r0 < a2.m; r0 += kPass2TileDim_) {
						const int mTile
						    = std::min(kPass2TileDim_, a2.m - r0);
						for (int c0 = 0; c0 < a2.n; c0 += kPass2TileDim_) {
							const int nTile
							    = std::min(kPass2TileDim_, a2.n - c0);

							long long colCursor = 0;
							for (SizeType c = 0; c < nChunks; ++c) {
								const int kc = std::max(
								    0,
								    std::min(chunkCols,
								             totalCols
								                 - static_cast<int>(
								                     colCursor)));
								// Clamp the column offset so
								// pointer arithmetic never strays
								// past this patch's allocated
								// column range (kc==0 chunks still
								// need a valid, in-bounds offset
								// even though no data will be
								// read/written there).
								const long long colOff = std::min(
								    colCursor,
								    static_cast<long long>(
								        totalCols));

								GemmArgs ac = a2;
								ac.m        = mTile;
								ac.n        = nTile;
								ac.k        = kc;
								ac.a_off = static_cast<long long>(
								               BXbatchOff[li])
								    + colOff
								        * static_cast<long long>(
								            ldB[ip])
								    + r0;
								ac.b_off = static_cast<long long>(
								               AbatchOff[li])
								    + colOff
								        * static_cast<long long>(
								            ldA[ip])
								    + c0;
								// ac.ldc keeps the FULL patch row
								// count (a2.ldc = rps[ip],
								// unchanged from above)
								// -- it is the real addressing
								// stride into this shard's d_vout,
								// distinct from ac.m (this tile's
								// row count) -- see matrixVector().
								// c_off is LOCAL to this shard's
								// d_vout (xyStart[ip] minus the
								// shard's own starting offset),
								// since unlike d_vin, d_vout is not
								// replicated.
								ac.c_off = static_cast<long long>(
								               xyStart[ip])
								    - static_cast<long long>(
								               xyStart[patchBegin])
								    + r0
								    + static_cast<long long>(c0)
								        * a2.ldc;
								pass2_args.push_back(ac);

								colCursor += chunkCols;
							}
						}
					}
				}

				shard.nbatch1 = pass1_args.size();
				shard.nbatch2 = pass2_args.size();
				shard.d_pass1 = DevArgsView(
				    Kokkos::view_alloc(shard.exec, "d_pass1"), shard.nbatch1);
				shard.d_pass2 = DevArgsView(
				    Kokkos::view_alloc(shard.exec, "d_pass2"), shard.nbatch2);
				{
					Kokkos::View<GemmArgs*,
					             Kokkos::HostSpace,
					             Kokkos::MemoryUnmanaged>
					    h1(pass1_args.data(), shard.nbatch1);
					Kokkos::View<GemmArgs*,
					             Kokkos::HostSpace,
					             Kokkos::MemoryUnmanaged>
					    h2(pass2_args.data(), shard.nbatch2);
					Kokkos::deep_copy(shard.exec, shard.d_pass1, h1);
					Kokkos::deep_copy(shard.exec, shard.d_pass2, h2);
				}

				shard.voutHostOffset = xyStart[patchBegin];
				shard.voutLocalSize  = xyStart[patchEnd] - xyStart[patchBegin];
				shard.d_vin          = ScalarView(
                                    Kokkos::view_alloc(
                                        shard.exec, Kokkos::WithoutInitializing, "d_vin"),
                                    xyStart[npatches]); // full replica of vin
				shard.d_vout = ScalarView(
				    Kokkos::view_alloc(
				        shard.exec, Kokkos::WithoutInitializing, "d_vout"),
				    shard.voutLocalSize);
			}

			// Execute ALL shards' pack tasks together in a single parallel
			// region on the host (packing is CPU-side and independent of
			// which GPU each task's data will eventually be copied to).
			if (!packTasks.empty()) {
				using HostExec = Kokkos::DefaultHostExecutionSpace;

				auto      tasksPtr = packTasks.data();
				const int ntasks   = static_cast<int>(packTasks.size());
				Kokkos::parallel_for(
				    "BatchedGemmKokkos::host_pack",
				    Kokkos::RangePolicy<HostExec>(0, ntasks),
				    KOKKOS_LAMBDA(const int idx) {
					    const PackTask& t = tasksPtr[idx];
					    // column-major copy: for each column j copy m entries
					    for (int j = 0; j < t.n; ++j) {
						    std::memcpy(
						        t.dst + static_cast<long long>(j) * t.ldd,
						        t.src + static_cast<long long>(j) * t.lds,
						        static_cast<size_t>(t.m)
						            * sizeof(ComplexOrRealType));
					    }
				    });
				Kokkos::fence();
			}
		}
		{
			Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup::rest");

			SizeType totalBatch1 = 0, totalBatch2 = 0;
			for (SizeType s = 0; s < numShards; ++s) {
				DeviceShard& shard = shards_[s];
				Kokkos::deep_copy(shard.exec, shard.d_flatAbatch, hostAbatch[s]);
				Kokkos::deep_copy(shard.exec, shard.d_flatBbatch, hostBbatch[s]);
				totalBatch1 += shard.nbatch1;
				totalBatch2 += shard.nbatch2;
			}
			for (SizeType s = 0; s < numShards; ++s)
				shards_[s].exec.fence();

			{
				PsimagLite::OstringStream msg(std::cout.precision());
				msg() << "setup done: npatches=" << npatches
				      << " noperator=" << noperator << " numDevices=" << numShards
				      << " pass1_batches=" << totalBatch1
				      << " pass2_batches=" << totalBatch2
				      << " totalAbatch=" << totalAbatch << "elems"
				      << " Bbatch=" << totalBbatch << "elems"
				      << " BXbatch=" << totalBXbatch << "elems";
				progress_.printline(msg, std::cout);
			}
		}
	}

private:

	// -----------------------------------------------------------------------

	const InitKronType&           initKron_;
	PsimagLite::ProgressIndicator progress_;
	mutable VectorMatrixType      garbage_; // owns sparse->dense expansions (freed in dtor)

	// One DeviceShard per physical device this instance's work is distributed
	// across.
	mutable std::vector<DeviceShard> shards_;
};

} // namespace Dmrg
#endif // BATCHED_GEMM_KOKKOS_H
