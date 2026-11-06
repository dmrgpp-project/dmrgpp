#ifndef ONE_TIME_EVOLUTION_HPP
#define ONE_TIME_EVOLUTION_HPP
#include <PsimagLite/AllocatorCpu.h>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Dmrg {

template <typename PvectorsType> class OneTimeEvolution {

public:

	using VectorWithOffsetType = typename PvectorsType::VectorWithOffsetType;
	using RealType             = typename PvectorsType::RealType;
	using VectorSizeType       = std::vector<SizeType>;
	using VectorRealType       = std::vector<RealType>;

	OneTimeEvolution(std::string srcKet,
	                 std::string name,
	                 SizeType    nsites,
	                 RealType    Eg,
	                 RealType    initialTime)
	    : timeParams_(parseTimeParams(std::move(name), nsites, Eg))
	    , indices_(timeParams_.timeSteps)
	    , srcKet_(std::move(srcKet))
	    , timesWithoutAdvancement_(0)
	    , time_(initialTime)
	{ }

	void
	initialize(SizeType firstIndex, const VectorWithOffsetType& src, PvectorsType& pVectors)
	{
		indices_[0] = firstIndex;
		for (SizeType i = 1; i < timeParams_.timeSteps; ++i) {
			auto lambda = [this, i](SizeType ind)
			{
				indices_[i] = ind;
				return std::string("|P") + std::to_string(ind) + ">";
			};

			pVectors.createNew(src, time_, lambda);
			if (i == timeParams_.timeSteps - 1)
				pVectors.registerLastKrylovSlot(firstIndex, indices_[i]);
		}
	}

	const VectorSizeType& indices() const { return indices_; }

	SizeType timesWithoutAdvancement() const { return timesWithoutAdvancement_; }

	RealType time() const { return time_; }

	SizeType timeSteps() const { return timeParams_.timeSteps; }

	RealType tau() const { return timeParams_.tau; }

	SizeType advanceEach() const { return timeParams_.advanceEach; }

	const std::string& algorithm() const { return timeParams_.algo; }

	const VectorRealType& chebyTransform() const { return timeParams_.chebyTransform; }

	SizeType disposition() const { return timeParams_.disposition; }

	const VectorSizeType& depends() const { return timeParams_.depends; }

	RealType energyForExp() const { return timeParams_.Eg; }

	void advanceTime(RealType tau) { time_ += tau; }

	void resetTimesWithoutAdvancement() { timesWithoutAdvancement_ = 1; }

	void incrementTimesWithoutAdvancement() { ++timesWithoutAdvancement_; }

	bool hasKet(SizeType ket) const
	{
		for (SizeType index : indices_)
			if (index == ket)
				return true;

		return false;
	}

	int sourceToDestroy() const
	{
		return (timeParams_.disposition != 2) ? -1 : getPindexOrMinusOne(srcKet_);
	}

private:

	struct TimeParams {

		TimeParams(SizeType nsites, RealType Eg1)
		    : timeSteps(3)
		    , tau(0.1)
		    , advanceEach(nsites - 2)
		    , algo("Krylov")
		    , disposition(0)
		    , Eg(Eg1)
		{ }

		SizeType       timeSteps;
		RealType       tau = 0.1;
		SizeType       advanceEach;
		std::string    algo;
		VectorRealType chebyTransform;
		SizeType       disposition;
		VectorSizeType depends;
		RealType       Eg;
	};

	static TimeParams parseTimeParams(std::string name, SizeType nsites, RealType Eg)
	{
		const std::string timeEvolve = "TimeEvolve";
		if (name.rfind(timeEvolve, 0) != 0)
			throw std::runtime_error("Only TimeEvolve implemented");

		TimeParams timeParams(nsites, Eg);
		name.erase(0, timeEvolve.length());
		std::string buffer;
		for (char c : name) {
			if (std::isspace(static_cast<unsigned char>(c)) || c == '{' || c == '}'
			    || c == '(' || c == ')')
				continue;
			buffer += c;
		}

		for (const std::string& token : split(buffer, ',')) {
			const auto equals = token.find('=');
			if (equals == std::string::npos
			    || token.find('=', equals + 1) != std::string::npos)
				throw std::runtime_error("Expected key=value in TimeEvolve");

			const std::string key   = token.substr(0, equals);
			const std::string value = token.substr(equals + 1);
			if (value == "~")
				continue;
			if (key == "advanceEach" || key == "AdvanceEach" || key == "advanceeach") {
				timeParams.advanceEach = std::stoul(value);
			} else if (key == "tau") {
				timeParams.tau = static_cast<RealType>(std::stod(value));
			} else if (key == "disposition") {
				timeParams.disposition = std::stoul(value);
			} else if (key == "steps") {
				timeParams.timeSteps = std::stoul(value);
			} else if (key == "algorithm") {
				timeParams.algo
				    = getChebyIfNeeded(timeParams.chebyTransform, value);
			} else if (key == "EnergyForExp" || key == "energyforexp") {
				timeParams.Eg = static_cast<RealType>(std::stod(value));
			} else if (key == "depends") {
				getDepends(timeParams.depends, value);
			} else {
				throw std::runtime_error("Unrecognized TimeEvolve key=" + key);
			}
		}

		return timeParams;
	}

	static std::vector<std::string> split(const std::string& input, char delimiter)
	{
		std::vector<std::string> result;
		SizeType                 begin = 0;
		while (true) {
			const SizeType end = input.find(delimiter, begin);
			result.push_back(input.substr(begin, end - begin));
			if (end == std::string::npos)
				break;
			begin = end + 1;
		}
		return result;
	}

	static std::string getChebyIfNeeded(VectorRealType& transform, const std::string& algo)
	{
		const std::string cheby = "Chebyshev";
		if (algo.rfind(cheby, 0) != 0)
			return algo;

		const std::vector<std::string> tokens = split(algo, ':');
		if (tokens.size() != 3 || tokens[0] != cheby)
			throw std::runtime_error("Expected Chebyshev:a:b");

		transform.resize(2);
		transform[0] = static_cast<RealType>(std::stod(tokens[1]));
		transform[1] = static_cast<RealType>(std::stod(tokens[2]));
		return algo;
	}

	static void getDepends(VectorSizeType& depends, const std::string& value)
	{
		for (const std::string& token : split(value, ':')) {
			if (token.length() < 2 || token[0] != 'P')
				throw std::runtime_error("Expected P followed by number");
			depends.push_back(std::stoul(token.substr(1)));
		}
	}

	static int getPindexOrMinusOne(std::string braOrKet)
	{
		if (braOrKet.length() > 2
		    && ((braOrKet.front() == '|' && braOrKet.back() == '>')
		        || (braOrKet.front() == '<' && braOrKet.back() == '|')))
			braOrKet = braOrKet.substr(1, braOrKet.length() - 2);

		if (braOrKet.empty() || (braOrKet[0] != 'P' && braOrKet[0] != 'R'))
			return -1;

		try {
			return static_cast<int>(std::stoul(braOrKet.substr(1)));
		} catch (const std::exception&) {
			return -1;
		}
	}

	TimeParams     timeParams_;
	VectorSizeType indices_;
	std::string    srcKet_;
	SizeType       timesWithoutAdvancement_;
	RealType       time_;
};
}
#endif // ONE_TIME_EVOLUTION_HPP