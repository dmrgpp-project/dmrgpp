#include "DensityOfStates.h"
#include "LatticeGf.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>

namespace {
using Complex = std::complex<double>;
using Dos     = Dmft::DensityOfStates<Complex>;
using Freq    = Dmft::FunctionOfFrequency<Complex>;

Complex semicircleGreen(Complex z, double halfBandwidth)
{
	const double w = z.imag();
	return Complex(0.0,
	               -2.0 * std::copysign(1.0, w)
	                   / (std::abs(w) + std::sqrt(w * w + halfBandwidth * halfBandwidth)));
}
}

TEST_CASE("semicircular DOS has the requested support, normalization, and moment",
          "[LatticeGf][DOS]")
{
	for (const double width : { 2.0, 4.0 }) {
		const double halfWidth = width / 2.0;
		const Dos    dos("semicircular", halfWidth);
		CHECK(dos.lowerBound() == -halfWidth);
		CHECK(dos.upperBound() == halfWidth);
		constexpr int n            = 200000;
		const double  dx           = width / n;
		double        norm         = 0.0;
		double        secondMoment = 0.0;
		for (int j = 0; j < n; ++j) {
			const double e   = -halfWidth + (j + 0.5) * dx;
			const double rho = dos(e);
			norm += rho * dx;
			secondMoment += e * e * rho * dx;
		}
		CHECK(norm == Catch::Approx(1.0).margin(2e-8));
		CHECK(secondMoment == Catch::Approx(width * width / 16.0).margin(2e-8));
	}
}

TEST_CASE("Bethe lattice GF obeys semicircle and self-consistency identities", "[LatticeGf]")
{
	constexpr double   beta  = 25.0;
	constexpr SizeType count = 80;
	for (const double width : { 2.0, 4.0 }) {
		const double t  = width / 4.0;
		const double t2 = t * t;
		Freq         sigma(beta, count);
		for (SizeType i = 0; i < sigma.totalMatsubaras(); ++i)
			sigma(i) = 0.0;
		Dmft::LatticeGf<Complex> lattice(
		    sigma, 0.0, "energy,semicircular," + std::to_string(width));
		lattice.update();
		for (SizeType i = 0; i < sigma.totalMatsubaras(); ++i) {
			const Complex z(0.0, sigma.omega(i));
			CHECK(lattice()(i).real()
			      == Catch::Approx(semicircleGreen(z, width / 2).real()).margin(1e-10));
			CHECK(lattice()(i).imag()
			      == Catch::Approx(semicircleGreen(z, width / 2).imag()).margin(1e-10));
		}

		for (SizeType i = 0; i < sigma.totalMatsubaras(); ++i) {
			const Complex iwn(0.0, sigma.omega(i));
			// A causal pole self-energy with the correct Matsubara conjugate symmetry.
			sigma(i) = 0.3 / (iwn - 0.2);
		}
		lattice.update();
		for (SizeType i = 0; i < sigma.totalMatsubaras(); ++i) {
			const Complex iwn(0.0, sigma.omega(i));
			const Complex g = lattice()(i);
			const Complex z = iwn - sigma(i);
			CHECK(std::abs(z * g - t2 * g * g - Complex(1.0, 0.0)) < 1e-8);
			CHECK(std::abs(1.0 / lattice.g0()(i) - (iwn - t2 * g)) < 1e-8);
		}
	}
}
