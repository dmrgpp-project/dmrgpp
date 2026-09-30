#ifndef FIT_H
#define FIT_H
#include "AndersonFunction.h"
#include "FitFunction.hpp"
#include "MinParams.h"
#include <PsimagLite/MersenneTwister.h>
#include <PsimagLite/Minimizer.h>
#include <PsimagLite/PsimagLite.h>

namespace Dmft {

template <typename ComplexOrRealType> class Fit {

public:

	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType          = std::vector<RealType>;
	using MinParamsType           = MinParams<RealType>;
	using AndersonFunctionType    = AndersonFunction<ComplexOrRealType>;
	using FitFunctionType         = FitFunction<ComplexOrRealType>;
	using FunctionOfFrequencyType = typename AndersonFunctionType::FunctionOfFrequencyType;
	using RngType                 = PsimagLite::MersenneTwister;
	using Options                 = typename FitFunctionType::Options;

	struct InitResults {

		InitResults(RealType ra_, RealType rb_, const VectorRealType& result_, bool reset_)
		    : ra(ra_)
		    , rb(rb_)
		    , result(result_)
		    , reset(reset_)
		{ }

		template <typename ReadableType>
		InitResults(ReadableType& io)
		    : ra(0)
		    , rb(0)
		    , reset(false)
		{
			try {
				io.readline(ra, "InitBathRa=");
			} catch (std::exception&) { }

			try {
				io.readline(rb, "InitBathRb=");
			} catch (std::exception&) { }

			try {
				io.read(result, "InitBathVector");
			} catch (std::exception&) { }

			try {
				int tmp = 0;
				io.readline(tmp, "InitBathReset=");
				reset = (tmp > 0);
			} catch (std::exception&) { }

			bool nonConstant = (result.size() > 0);
			bool isConstant  = (ra != 0 || rb != 0);
			if (isConstant && nonConstant)
				err("InitBath: Cannot have ra or rb and also a vector of init "
				    "results\n");

			if (!nonConstant && !isConstant)
				ra = 1;
		}

		RealType       ra;
		RealType       rb;
		VectorRealType result;
		bool           reset;
	};

	Fit(SizeType nBath, const MinParamsType& minParams, const InitResults& initResults)
	    : nBath_(nBath)
	    , minParams_(minParams)
	    , results_(2 * nBath_)
	    , rng_(1234)
	    , initResults_(initResults)
	{ }

	// Compute the optimized bath parameters and store them in vector g0
	// See AndersonFunction.h documentation for the fit function, and
	// for the order of storage of bath parameters
	void fit(const FunctionOfFrequencyType& g0, const RealType& mu, Options options)
	{
		FitFunctionType f(nBath_, g0, mu, options);

		VectorRealType results(f.size());
		setResults(results);

		using MinimizerType = PsimagLite::Minimizer<RealType, FitFunctionType>;
		MinimizerType min(f, minParams_.maxIter, minParams_.verbose);
		int iter = 0;
		if (minParams_.method == MinParamsType::Method::CONJUGATE_GRADIENT) {
			iter = min.conjugateGradient(
			    results, minParams_.delta, minParams_.delta2, minParams_.tolerance);

			if (minimizerFailed(min)) {
				if (PsimagLite::Concurrency::root()) {
					std::cerr << "Bath fit conjugate-gradient stalled after "
					          << completedIterations(iter) << " iterations: GSL status "
					          << min.status() << " (" << gsl_strerror(min.status())
					          << "); retrying simplex from the warm bath\n";
				}

				iter = min.simplex(results, minParams_.delta, minParams_.tolerance);
				throwIfMinimizerFailed(min, "simplex", iter);
			}
		} else {
			iter = min.simplex(results, minParams_.delta, minParams_.tolerance);
			throwIfMinimizerFailed(min, "simplex", iter);
		}

		assert(results.size() == nBath_ || results.size() == 2 * nBath_);
		if (results.size() == 2 * nBath_) {
			for (SizeType i = 0; i < 2 * nBath_; ++i) {
				results_[i] = results[i];
			}
		} else {
			// particle-hole symmetric case

			// The energy zero bath site is the middle one (nBath odd case only)

			// copy Vs first, starting with the fitted ones...
			SizeType number_of_fitted_Vs = (nBath_ & 1) ? (nBath_ + 1) / 2 : nBath_ / 2;
			for (SizeType i = 0; i < number_of_fitted_Vs; ++i) {
				results_[i] = results[i];
			}

			// ...and then the rest of the Vs are a mirror
			// Works also for Nbath_ odd
			for (SizeType i = number_of_fitted_Vs; i < nBath_; ++i) {
				results_[i] = results[i - number_of_fitted_Vs];
			}

			// copy onsite energies, first the fitted ones...
			SizeType offset1 = nBath_ - number_of_fitted_Vs;
			for (SizeType i = number_of_fitted_Vs; i < nBath_; ++i) {
				results_[i + offset1] = results[i];
			}

			// .. then the center bath site (if nbath is odd)...
			SizeType one_or_zero = (nBath_ & 1);
			if (nBath_ & 1) {
				results_[nBath_ + offset1] = 0;
			}

			// ...and the rest are the opposites
			// Works also for Nbath_ odd
			SizeType offset2 = 2 * (nBath_ - number_of_fitted_Vs);
			for (SizeType i = number_of_fitted_Vs; i < nBath_; ++i) {
				results_[i + offset2 + one_or_zero] = -results[i];
			}

			assert(nBath_ + offset2 + one_or_zero == 2 * nBath_);
		}

		hasFitted_ = true;
	}

	const VectorRealType& result() const { return results_; }

	SizeType nBath() const { return nBath_; }

	static Options computeOptions(const std::string& options)
	{
		return FitFunctionType::computeOptions(options);
	}

private:

	template <typename MinimizerType>
	bool minimizerFailed(const MinimizerType& min) const
	{
		return minParams_.maxIter > 0 && min.status() != MinimizerType::GSL_SUCCESS
		       && min.status() != MinimizerType::GSL_CONTINUE;
	}

	static int completedIterations(int iter) { return (iter < 0) ? -iter : iter; }

	template <typename MinimizerType>
	void throwIfMinimizerFailed(const MinimizerType&       min,
	                            const PsimagLite::String& method,
	                            int                         iter) const
	{
		if (!minimizerFailed(min))
			return;

		PsimagLite::String msg("Bath fit failed with ");
		msg += method + " after " + ttos(completedIterations(iter)) + " iterations: GSL status "
		    + ttos(min.status()) + " (" + gsl_strerror(min.status()) + ")\n";
		err(msg);
	}

	void setResults(VectorRealType& results)
	{
		if (hasFitted_ && !initResults_.reset) {
			setResultsFromPreviousFit(results);
			return;
		}

		bool nonConstant = (initResults_.result.size() > 0);
		bool isConstant  = (initResults_.ra != 0 || initResults_.rb != 0);
		if (isConstant && nonConstant)
			err("InitResults: Cannot have ra or rb and also a vector of init "
			    "results\n");

		if (nonConstant && initResults_.result.size() != results.size())
			err(PsimagLite::String(
			        "InitResults: vector of init results has wrong size: ")
			    + "expected " + ttos(results.size()) + ", but found "
			    + ttos(initResults_.result.size()) + "\n");

		for (SizeType i = 0; i < results.size(); ++i)
			results[i] = (isConstant) ? initResults_.ra * rng_() + initResults_.rb
			                          : initResults_.result[i];
	}

	void setResultsFromPreviousFit(VectorRealType& results) const
	{
		assert(results_.size() == 2 * nBath_);
		if (results.size() == 2 * nBath_) {
			results = results_;
			return;
		}

		assert(results.size() == nBath_);
		SizeType numberOfIndependentVs = (nBath_ & 1) ? (nBath_ + 1) / 2 : nBath_ / 2;
		for (SizeType i = 0; i < numberOfIndependentVs; ++i)
			results[i] = results_[i];

		SizeType energyOffset = nBath_ - numberOfIndependentVs;
		for (SizeType i = numberOfIndependentVs; i < nBath_; ++i)
			results[i] = results_[i + energyOffset];
	}

	const SizeType       nBath_; // number of bath sites
	const MinParamsType& minParams_; // parameters for fitting algorithm
	VectorRealType       results_; // stores bath parameters
	bool                 hasFitted_ = false;
	RngType              rng_;
	const InitResults&   initResults_;
};
}
#endif // FIT_H
