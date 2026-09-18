// Standalone, deterministic microbenchmark for BatchedGemm::matrixVector.
// Unlike running the full dmrg solver on a real input, this bypasses the Lanczos eigensolver
// entirely -- so the exact same fixed synthetic Kron problem (patch count/sizes, connection count)
// is exercised on every call, with no run-to-run variation in which quantum number sectors get
// explored.
//
// Usage: bench_BatchedGemm [npatches] [nSmall] [smallSize] [nLarge] [largeSize] [iters]
//   Builds `npatches` independent (block-diagonal) patches: `nSmall` of them
//   sized `smallSize` x `smallSize` and the rest ((npatches - nSmall), or
//   `nLarge` if given explicitly) sized `largeSize` x `largeSize` -- modeling
//   the severe task-size imbalance (a few huge patches among many small
//   ones) observed in the real 345 workload's setup_() diagnostics.

#include "DMRGConfig.h"
#include "MatrixVectorKron/BatchedGemmInclude.hh"
#include <Kokkos_Core.hpp>
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Vector.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using SizeType = unsigned long;

// --- Minimal fake types (same shape as dmrg/Engine/tests/test_BatchedGemm.cpp) ---

struct FakeParams {
	struct Opt {
		bool isSet(const std::string& s) const { return s == "BatchedGemm"; }
	} options;
};

struct FakeLR {
	std::vector<int> parts;
	int              partition(size_t i) const { return parts[i]; }
	int              size() const { return parts.back(); }
};

struct FakeLrs {
	FakeLR        l, r;
	const FakeLR& left() const { return l; }
	const FakeLR& right() const { return r; }
};

struct FakeMatrixDenseOrSparse {
	using value_type = double;
	using VectorType = PsimagLite::Vector<value_type>::Type;
	using MatrixType = PsimagLite::Matrix<value_type>;
	using CrsType    = PsimagLite::CrsMatrix<value_type>;
	MatrixType mat;
	CrsType    crs;
	bool       zero = false;
	FakeMatrixDenseOrSparse()
	    : mat()
	    , crs(0, 0)
	    , zero(true)
	{ }
	FakeMatrixDenseOrSparse(const MatrixType& m)
	    : mat(m)
	    , crs(0, 0)
	    , zero(false)
	{ }
	bool              isZero() const { return zero; }
	bool              isDense() const { return true; }
	const MatrixType& dense() const { return mat; }
	const CrsType&    sparse() const { return crs; }
	int               rows() const { return mat.rows(); }
	int               cols() const { return mat.cols(); }
};

struct FakeArrayOfMatStructType {
	using MatrixDenseOrSparseType = FakeMatrixDenseOrSparse;
	std::vector<std::vector<FakeMatrixDenseOrSparse>> storage;
	FakeArrayOfMatStructType(size_t npatches = 0)
	{
		storage.resize(npatches);
		for (size_t i = 0; i < npatches; ++i)
			storage[i].resize(npatches);
	}
	const FakeMatrixDenseOrSparse* operator()(size_t ip, size_t jp) const
	{
		if (storage[ip][jp].isZero())
			return nullptr;
		return &storage[ip][jp];
	}
	FakeMatrixDenseOrSparse* operator()(size_t ip, size_t jp)
	{
		if (storage[ip][jp].isZero())
			return nullptr;
		return &storage[ip][jp];
	}
	void set(size_t ip, size_t jp, const FakeMatrixDenseOrSparse& m) { storage[ip][jp] = m; }
};

struct GenIjPatch {
	enum LeftOrRightEnumType
	{
		LEFT  = 0,
		RIGHT = 1
	};
	using BasisType = int;
};

struct FakeInitKron {
	enum WhatBasisEnum
	{
		OLD = 0,
		NEW = 1
	};

	using ArrayOfMatStructType = FakeArrayOfMatStructType;
	using GenIjPatchType       = GenIjPatch;
	using SparseMatrixType     = PsimagLite::CrsMatrix<double>;

	FakeParams               p;
	FakeArrayOfMatStructType xc0, yc0;
	FakeLrs                  lrs_;
	std::vector<size_t>      patchLeft, patchRight;
	size_t                   npatches_  = 0;
	size_t                   noperator_ = 0;

	FakeInitKron(size_t npatches, size_t noperator)
	    : p()
	    , xc0(npatches)
	    , yc0(npatches)
	    , npatches_(npatches)
	    , noperator_(noperator)
	{
		patchLeft.resize(npatches_, 0);
		patchRight.resize(npatches_, 0);
	}

	const FakeParams& params() const { return p; }
	size_t            numberOfPatches(int) const { return npatches_; }
	size_t            connections() const { return noperator_; }

	const std::vector<size_t>& patch(WhatBasisEnum, GenIjPatch::LeftOrRightEnumType side) const
	{
		if (side == GenIjPatch::LEFT)
			return patchLeft;
		return patchRight;
	}

	const FakeLrs& lrs(WhatBasisEnum) const { return lrs_; }

	SizeType offsetForPatches(WhatBasisEnum, size_t ipatch) const
	{
		SizeType off = 0;
		for (size_t p = 0; p < ipatch; ++p) {
			size_t igroup = patchLeft[p];
			size_t jgroup = patchRight[p];
			int    L1     = lrs_.left().partition(igroup);
			int    L2     = lrs_.left().partition(igroup + 1);
			int    R1     = lrs_.right().partition(jgroup);
			int    R2     = lrs_.right().partition(jgroup + 1);
			off += static_cast<SizeType>((L2 - L1) * (R2 - R1));
		}
		return off;
	}

	const FakeArrayOfMatStructType& xc(size_t) const { return xc0; }
	const FakeArrayOfMatStructType& yc(size_t) const { return yc0; }

	void
	checks(const FakeMatrixDenseOrSparse&, const FakeMatrixDenseOrSparse&, size_t, size_t) const
	{ }
};

namespace {

PsimagLite::Matrix<double> randomSquare(int n, std::mt19937_64& rng)
{
	std::uniform_real_distribution<double> dist(-1.0, 1.0);
	PsimagLite::Matrix<double>             m(n, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			m(i, j) = dist(rng);
	return m;
}

} // namespace

int main(int argc, char** argv)
{
	Kokkos::ScopeGuard scope_guard(argc, argv);

	// Defaults chosen to mimic the severe task-size imbalance observed in
	// setup_()'s host_pack diagnostics for input345.ain: mostly small
	// patches with a handful of much larger ones.
	const int npatches  = (argc > 1) ? std::atoi(argv[1]) : 40;
	const int nSmall    = (argc > 2) ? std::atoi(argv[2]) : 34;
	const int smallSize = (argc > 3) ? std::atoi(argv[3]) : 60;
	const int nLarge    = (argc > 4) ? std::atoi(argv[4]) : (npatches - nSmall);
	const int largeSize = (argc > 5) ? std::atoi(argv[5]) : 900;
	const int iters     = (argc > 6) ? std::atoi(argv[6]) : 30;

	if (nSmall + nLarge != npatches) {
		std::cerr << "nSmall + nLarge must equal npatches\n";
		return 1;
	}

	std::mt19937_64 rng(42);

	FakeInitKron fk(static_cast<size_t>(npatches), 1);

	std::vector<int> sizes(npatches);
	for (int i = 0; i < npatches; ++i)
		sizes[i] = (i < nSmall) ? smallSize : largeSize;
	// Shuffle so large/small patches are interleaved rather than grouped
	// (grouping would artificially favor static-scheduling load balance).
	std::shuffle(sizes.begin(), sizes.end(), rng);

	fk.lrs_.l.parts.resize(npatches + 1);
	fk.lrs_.r.parts.resize(npatches + 1);
	fk.lrs_.l.parts[0] = 0;
	fk.lrs_.r.parts[0] = 0;
	for (int i = 0; i < npatches; ++i) {
		fk.lrs_.l.parts[i + 1] = fk.lrs_.l.parts[i] + sizes[i];
		fk.lrs_.r.parts[i + 1] = fk.lrs_.r.parts[i] + sizes[i];
		fk.patchLeft[i]        = static_cast<size_t>(i);
		fk.patchRight[i]       = static_cast<size_t>(i);
		fk.xc0.set(i, i, FakeMatrixDenseOrSparse(randomSquare(sizes[i], rng)));
		fk.yc0.set(i, i, FakeMatrixDenseOrSparse(randomSquare(sizes[i], rng)));
	}

	SizeType totalXY = 0;
	for (int i = 0; i < npatches; ++i)
		totalXY += static_cast<SizeType>(sizes[i]) * static_cast<SizeType>(sizes[i]);

	std::cout << "Synthetic problem: npatches=" << npatches << " (small=" << nSmall << "x"
	          << smallSize << ", large=" << nLarge << "x" << largeSize
	          << "), total vin/vout size=" << totalXY << "\n";

	using VectorType = PsimagLite::Vector<double>::Type;
	VectorType vin(totalXY), vout(totalXY);
	{
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		for (SizeType i = 0; i < totalXY; ++i)
			vin[i] = dist(rng);
	}

	auto                             t0 = std::chrono::steady_clock::now();
	Dmrg::BATCHED_GEMM<FakeInitKron> bgk(fk);
	auto                             t1 = std::chrono::steady_clock::now();
	std::cout << "setup_() time: " << std::chrono::duration<double, std::milli>(t1 - t0).count()
	          << " ms\n";

	// Warm-up call (not timed): first-touch / lazy allocation effects, JIT,
	// etc. should not pollute the measured iterations below.
	bgk.matrixVector(vout, vin);

	std::vector<double> callTimesMs;
	callTimesMs.reserve(iters);
	for (int it = 0; it < iters; ++it) {
		auto s = std::chrono::steady_clock::now();
		bgk.matrixVector(vout, vin);
		auto e = std::chrono::steady_clock::now();
		callTimesMs.push_back(std::chrono::duration<double, std::milli>(e - s).count());
	}

	std::sort(callTimesMs.begin(), callTimesMs.end());
	const double sum = std::accumulate(callTimesMs.begin(), callTimesMs.end(), 0.0);
	const double avg = sum / callTimesMs.size();
	const double med = callTimesMs[callTimesMs.size() / 2];
	std::cout << "matrixVector: iters=" << iters << " min=" << callTimesMs.front()
	          << " ms median=" << med << " ms avg=" << avg << " ms max=" << callTimesMs.back()
	          << " ms\n";

	return 0;
}
