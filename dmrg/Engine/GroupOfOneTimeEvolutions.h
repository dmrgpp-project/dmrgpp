#ifndef GROUP_OF_ONE_TIME_EVOLUTIONS_H
#define GROUP_OF_ONE_TIME_EVOLUTIONS_H
#include "OneTimeEvolution.hpp"
#include <vector>

namespace Dmrg {

template <typename PvectorsType> class GroupOfOneTimeEvolutions {

public:

	using OneTimeEvolutionType       = OneTimeEvolution<PvectorsType>;
	using VectorOneTimeEvolutionType = std::vector<OneTimeEvolutionType*>;
	using RealType                   = typename PvectorsType::RealType;

	GroupOfOneTimeEvolutions() { }

	GroupOfOneTimeEvolutions(const GroupOfOneTimeEvolutions&) = delete;

	GroupOfOneTimeEvolutions& operator=(const GroupOfOneTimeEvolutions&) = delete;

	~GroupOfOneTimeEvolutions()
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete vEvolutions_[i];
			vEvolutions_[i] = nullptr;
		}
	}

	OneTimeEvolutionType* findThisEvolution(SizeType firstIndex) const
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i)
			if (vEvolutions_[i]->indices()[0] == firstIndex)
				return vEvolutions_[i];

		return nullptr;
	}

	void pushBack(OneTimeEvolutionType* ptr) { vEvolutions_.push_back(ptr); }

	int getLastIndex(SizeType firstIndex) const
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i) {
			const auto& idx = vEvolutions_[i]->indices();
			if (idx[0] == firstIndex)
				return idx[idx.size() - 1];
		}
		return -1;
	}

	SizeType size() const { return vEvolutions_.size(); }

	RealType getTimeForKet(SizeType ket) const
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i) {
			if (vEvolutions_[i]->hasKet(ket))
				return vEvolutions_[i]->time();
		}

		return 0;
	}

private:

	VectorOneTimeEvolutionType vEvolutions_;
};
}
#endif // GROUP_OF_ONE_TIME_EVOLUTIONS_H