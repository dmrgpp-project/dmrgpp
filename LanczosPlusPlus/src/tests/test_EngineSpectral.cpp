#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/Engine.h"
#include "LanczosPlusPlus/src/Engine/InputCheck.h"
#include "LanczosPlusPlus/src/Engine/InternalProductOnTheFly.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include <PsimagLite/Geometry/Geometry.h>
#include <PsimagLite/InputNg.h>
#include <PsimagLite/Matsubaras.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/RealFrequencyRange.hh>
#include <array>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>

int main(int argc, char** argv)
{
	PsimagLite::PsiApp app("test_EngineSpectral", &argc, &argv, 1);
	return Catch::Session().run(argc, argv);
}

namespace {
using Complex              = std::complex<double>;
using Collection           = PsimagLite::ContinuedFractionCollection<double>;
using Operator             = LanczosPlusPlus::LabeledOperator;
using Input                = PsimagLite::InputNg<LanczosPlusPlus::InputCheck>;
constexpr double tolerance = 2e-8;
constexpr double hopping   = 0.7;
constexpr double e0        = -0.3;
constexpr double e1        = 0.4;

const PsimagLite::String hubbardInput = R"(
TotalNumberOfSites=2
NumberOfTerms=1
DegreesOfFreedom=1
GeometryKind=chain
GeometryOptions=ConstantValues
Connectors 1 0.7
Model=HubbardOneBand
hubbardU 2 0 0
potentialV 4 -0.3 0.4 -0.3 0.4
SolverOptions=none
TargetElectronsUp=1
TargetElectronsDown=1
)";

const PsimagLite::String twoOrbitalInput = R"(
TotalNumberOfSites=2
NumberOfTerms=1
DegreesOfFreedom=2
GeometryKind=chain
GeometryOptions=ConstantValues
Connectors 2 2
0.7 0.2
0.2 0.5
Model=FeAsBasedSc
Orbitals=2
FeAsMode=INT_PAPER33
hubbardU 4 0 0 0 0
potentialV 8 -0.3 0.4 0.2 -0.1 -0.3 0.4 0.2 -0.1
SolverOptions=none
TargetElectronsUp=1
TargetElectronsDown=1
)";

// The same concrete model/engine stack as ImpuritySolverExactDiag, without
// any DMFT normalization, tail fitting, or postprocessing.
template <typename Scalar> struct System {
	using Geometry
	    = PsimagLite::Geometry<Scalar, Input::Readable, LanczosPlusPlus::LanczosGlobals>;
	using Selector = LanczosPlusPlus::ModelSelector<Scalar, Geometry, Input::Readable>;
	using Model    = typename Selector::ModelBaseType;
	using Symmetry = LanczosPlusPlus::DefaultSymmetry<Geometry, typename Model::BasisBaseType>;
	using InternalProduct = LanczosPlusPlus::InternalProductOnTheFly<Model, Symmetry>;
	using Engine
	    = LanczosPlusPlus::Engine<Model, LanczosPlusPlus::InternalProductOnTheFly, Symmetry>;

	explicit System(const PsimagLite::String& text)
	    : writable(check, text)
	    , readable(writable)
	    , geometry(readable)
	    , selector(readable, geometry)
	    , engine(selector(), readable)
	{ }

	LanczosPlusPlus::InputCheck check;
	Input::Writeable            writable;
	Input::Readable             readable;
	Geometry                    geometry;
	Selector                    selector;
	Engine                      engine;
};

Complex frequency(const PsimagLite::RealFrequencyRange<double>& grid, SizeType i)
{
	return { grid.omega(i), grid.delta() };
}

Complex frequency(const PsimagLite::Matsubaras<double>& grid, SizeType i)
{
	return { grid.delta(), grid.omega(i) };
}

template <typename Function> void bothAxes(Function test)
{
	SECTION("REAL")
	{
		test(PsimagLite::FreqEnum::REAL,
		     PsimagLite::RealFrequencyRange<double>(-2.1, 0.37, 13, 0.2));
	}
	SECTION("MATSUBARA")
	{
		test(PsimagLite::FreqEnum::MATSUBARA, PsimagLite::Matsubaras<double>(7., 6, 0.));
	}
}

void checkNear(Complex actual, Complex expected)
{
	INFO("actual=" << actual << " expected=" << expected);
	CHECK(std::abs(actual - expected) < tolerance);
}

template <typename Scalar>
Collection spectrum(System<Scalar>&      system,
                    PsimagLite::FreqEnum axis,
                    Operator::Label      label,
                    SizeType             i,
                    SizeType             j,
                    SizeType             orbI = 0,
                    SizeType             orbJ = 0)
{
	Collection                                        result(axis);
	typename System<Scalar>::Engine::VectorStringType labels;
	system.engine.spectralFunction(result,
	                               labels,
	                               Operator(label),
	                               i,
	                               j,
	                               std::make_pair(0, 0),
	                               std::make_pair(orbI, orbJ));
	CHECK(labels.size() == result.size());
	return result;
}

// Direct inverse of zI-h; no Engine Hamiltonian, eigenvectors, or CF data
// enter the reference calculation. Partial-pivoted Gauss-Jordan suffices
// for these tiny matrices at frequencies strictly away from the poles.
template <SizeType N>
std::array<std::array<Complex, N>, N> resolvent(const std::array<std::array<double, N>, N>& h,
                                                Complex                                     z)
{
	std::array<std::array<Complex, 2 * N>, N> a {};
	for (SizeType i = 0; i < N; ++i) {
		for (SizeType j = 0; j < N; ++j)
			a[i][j] = (i == j ? z : Complex(0)) - h[i][j];
		a[i][N + i] = 1.;
	}
	for (SizeType col = 0; col < N; ++col) {
		SizeType pivot = col;
		for (SizeType row = col + 1; row < N; ++row)
			if (std::abs(a[row][col]) > std::abs(a[pivot][col]))
				pivot = row;
		std::swap(a[col], a[pivot]);
		REQUIRE(std::abs(a[col][col]) > 1e-12);
		const Complex divisor = a[col][col];
		for (auto& value : a[col])
			value /= divisor;
		for (SizeType row = 0; row < N; ++row) {
			if (row == col)
				continue;
			const Complex factor = a[row][col];
			for (SizeType j = 0; j < 2 * N; ++j)
				a[row][j] -= factor * a[col][j];
		}
	}
	std::array<std::array<Complex, N>, N> inverse {};
	for (SizeType i = 0; i < N; ++i)
		for (SizeType j = 0; j < N; ++j)
			inverse[i][j] = a[i][N + j];
	return inverse;
}

const std::array<std::array<double, 2>, 2> dimerH { { { e0, hopping }, { hopping, e1 } } };
const double gap   = std::sqrt((e0 - e1) * (e0 - e1) + 4 * hopping * hopping);
const double lower = (e0 + e1 - gap) / 2;
const double upper = (e0 + e1 + gap) / 2;

template <typename Scalar> void checkFermion(bool local)
{
	System<Scalar> system(hubbardInput);
	CHECK(std::abs(system.engine.energies(0) - 2 * lower) < tolerance);
	bothAxes(
	    [&](auto axis, const auto& grid)
	    {
		    const SizeType count = local ? 2 : 1;
		    for (SizeType site = 0; site < count; ++site) {
			    CAPTURE(site);
			    auto cf = spectrum(
			        system, axis, Operator::Label::OPERATOR_C, site, local ? site : 1);
			    REQUIRE(cf.size() == (local ? 2 : 4));
			    Collection::PlotDataType total;
			    cf.plot(total, grid);
			    REQUIRE(total.size() == grid.total());
			    Collection::PlotDataType particle, hole;
			    if (local) {
				    cf.plotOne(0, particle, grid);
				    cf.plotOne(1, hole, grid);
				    REQUIRE(particle.size() == grid.total());
				    REQUIRE(hole.size() == grid.total());
			    }
			    for (SizeType k = 0; k < grid.total(); ++k) {
				    const auto z = frequency(grid, k);
				    CAPTURE(z);
				    CHECK(total[k].first == grid.omega(k));
				    const auto g = resolvent(dimerH, z);
				    // Four-branch polarization returns G01+G10, not G01.
				    checkNear(total[k].second,
				              local ? g[site][site] : g[0][1] + g[1][0]);
				    if (local) {
					    const double occupied
					        = (upper - (site == 0 ? e0 : e1)) / gap;
					    // Unit total spectral weight: these residues add to
					    // one.
					    checkNear(particle[k].second,
					              (1 - occupied) / (z - upper));
					    checkNear(hole[k].second, occupied / (z - lower));
				    }
			    }
		    }
	    });
}

template <typename Scalar> void checkInterorbital()
{
	System<Scalar> system(twoOrbitalInput);
	// FeBasedSc potential index: site + (orbital + orbitals*spin)*nsite.
	// Reference mode order here is (site0/orb0, site0/orb1,
	// site1/orb0, site1/orb1), NOT the potential input's orbital-major order.
	const std::array<std::array<double, 4>, 4> h { { { -0.3, 0., 0.7, 0.2 },
		                                         { 0., 0.2, 0.2, 0.5 },
		                                         { 0.7, 0.2, 0.4, 0. },
		                                         { 0.2, 0.5, 0., -0.1 } } };
	bothAxes(
	    [&](auto axis, const auto& grid)
	    {
		    auto cf = spectrum(system, axis, Operator::Label::OPERATOR_C, 0, 0, 0, 1);
		    REQUIRE(cf.size() == 4);
		    Collection::PlotDataType total;
		    cf.plot(total, grid);
		    REQUIRE(total.size() == grid.total());
		    for (SizeType k = 0; k < grid.total(); ++k) {
			    const auto z = frequency(grid, k);
			    CAPTURE(z);
			    const auto g = resolvent(h, z);
			    checkNear(total[k].second, g[0][1] + g[1][0]);
		    }
	    });
}

template <typename Scalar> void checkDensity(bool local)
{
	System<Scalar> system(hubbardInput);
	CHECK(std::abs(system.engine.energies(0) - 2 * lower) < tolerance);
	bothAxes(
	    [&](auto axis, const auto& grid)
	    {
		    // OPERATOR_N with spin=0 is n_up, not n_up+n_down.
		    auto cf = spectrum(system, axis, Operator::Label::OPERATOR_N, 0, local ? 0 : 1);
		    REQUIRE(cf.size() == (local ? 2 : 4));
		    Collection::PlotDataType total;
		    cf.plot(total, grid);
		    REQUIRE(total.size() == grid.total());
		    const double weight = hopping * hopping / (gap * gap);
		    for (SizeType k = 0; k < grid.total(); ++k) {
			    const auto z = frequency(grid, k);
			    CAPTURE(z);
			    const auto chi = weight * (1. / (z - gap) - 1. / (z + gap));
			    checkNear(total[k].second, local ? chi : -2. * chi);
		    }
	    });
}
} // namespace

TEST_CASE("Engine raw local fermion has unit spectral weight", "[EngineSpectral][local]")
{
	SECTION("real scalar") { checkFermion<double>(true); }
	SECTION("complex scalar") { checkFermion<Complex>(true); }
}

TEST_CASE("Engine raw offsite fermion is G01 plus G10", "[EngineSpectral][offsite]")
{
	SECTION("real scalar") { checkFermion<double>(false); }
	SECTION("complex scalar") { checkFermion<Complex>(false); }
}

TEST_CASE("Engine same-site different-orbital fermion retains four branches",
          "[EngineSpectral][interorbital]")
{
	SECTION("real scalar") { checkInterorbital<double>(); }
	SECTION("complex scalar") { checkInterorbital<Complex>(); }
}

TEST_CASE("Engine raw local density response", "[EngineSpectral][local][density]")
{
	SECTION("real scalar") { checkDensity<double>(true); }
	SECTION("complex scalar") { checkDensity<Complex>(true); }
}

TEST_CASE("Engine raw offsite density is the symmetric response",
          "[EngineSpectral][offsite][density]")
{
	SECTION("real scalar") { checkDensity<double>(false); }
	SECTION("complex scalar") { checkDensity<Complex>(false); }
}
