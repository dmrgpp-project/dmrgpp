// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup PsimagLite */
/*@{*/

/*! \file LadderBath.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LADDER_BATH_H
#define LADDER_BATH_H

#include "Ladder.h"

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType>
class LadderBath : public GeometryBase<ComplexOrRealType, InputType> {

	using PairType   = std::pair<int, int>;
	using LadderType = Ladder<ComplexOrRealType, InputType>;

public:

	enum
	{
		DIRECTION_X = LadderType::DIRECTION_X,
		DIRECTION_Y = LadderType::DIRECTION_Y,
		DIRECTION_BATH
	};

	LadderBath() { }

	LadderBath(SizeType linSize, InputType& io)
	    : linSize_(linSize)
	    , ladder_(0)
	{
		io.readline(bathSitesPerSite_, "BathSitesPerSite=");
		clusterSize_ = linSize_ / (1 + bathSitesPerSite_);

		ladder_ = new LadderType(clusterSize_, io);
	}

	~LadderBath() override
	{
		if (ladder_)
			delete ladder_;
	}

	SizeType dirs() const override { return 3; }

	SizeType length(SizeType) const override { return this->unimplemented("length"); }

	SizeType translate(SizeType, SizeType, SizeType) const override
	{
		return this->unimplemented("translate");
	}

	SizeType getVectorSize(SizeType dirId) const override
	{
		if (dirId == DIRECTION_BATH)
			return bathSitesPerSite_ * clusterSize_;
		return ladder_->getVectorSize(dirId);
	}

	bool connected(SizeType i1, SizeType i2) const override
	{
		if (i1 == i2)
			return false;
		int c1 = getClusterSite(i1).first;
		int c2 = getClusterSite(i2).first;
		// 4 possibilites
		//  1. both in the cluster
		if (c1 < 0 && c2 < 0)
			return connectedInCluster(i1, i2);
		// both in the bath:
		if (c1 >= 0 && c2 >= 0)
			return false;
		// cluster - bath
		if (c1 < 0) {
			return (SizeType(c2) == i1) ? true : false;
		}
		// bath - cluster
		return (SizeType(c1) == i2) ? true : false;
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1, SizeType i2) const override
	{
		int c1 = getClusterSite(i1).first;
		int c2 = getClusterSite(i2).first;
		// two possibilities
		// 1. both in the cluster
		if (c1 < 0 && c2 < 0)
			return calcDirInCluster(i1, i2);
		// cluster - bath or bath cluster:
		return DIRECTION_BATH;
	}

	bool fringe(SizeType i, SizeType smax, SizeType emin) const override
	{
		int c = getClusterSite(i).first;
		if (c >= 0)
			return false; // no bath site is ever fringe
		return fringeInCluster(i, smax, emin);
	}

	// assumes i1 and i2 are connected
	SizeType handle(SizeType i1, SizeType i2) const override
	{
		PairType c1 = getClusterSite(i1);
		PairType c2 = getClusterSite(i2);
		// two possibilities
		// 1. both in the cluster
		if (c1.first < 0 && c2.first < 0)
			return handleInCluster(i1, i2);
		// cluster - bath or bath cluster
		PairType x = (c1.first < 0) ? c2 : c1;
		if (x.first < 0 || x.second < 0)
			throw RuntimeError("Internal error in handle\n");
		SizeType firstClusterSite = (clusterSize_ / 2) * bathSitesPerSite_;
		x.first -= firstClusterSite;

		return x.first * bathSitesPerSite_ + x.second;
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax, SizeType emin, SizeType siteNew) const override
	{
		PairType c1 = getClusterSite(siteNew);

		// in the cluster
		if (c1.first < 0) {
			SizeType firstClusterSite = (clusterSize_ / 2) * bathSitesPerSite_;
			SizeType siteNewCluster   = siteNew - firstClusterSite;
			SizeType smaxCluster      = smax - firstClusterSite;
			SizeType eminCluster      = emin - firstClusterSite;
			SizeType siteSubs
			    = ladder_->getSubstituteSite(smaxCluster, eminCluster, siteNewCluster);
			return siteSubs + firstClusterSite;
		}

		String str(__FILE__);
		str += " " + ttos(__LINE__) + "Internal error in getSubstituteSite\n";
		throw RuntimeError(str);
	}

	String label() const override { return "ladderbath"; }

	SizeType maxConnections() const override { return clusterSize_ + 1; }

	SizeType findReflection(SizeType) const override
	{
		throw RuntimeError("findReflection: unimplemented (sorry)\n");
	}

private:

	// if i is in the cluster return the pair (-1,-1)
	// else return the corresponding cluster site c and the number
	// of this bath site as a pair (c,b)
	PairType getClusterSite(SizeType i) const
	{
		SizeType firstClusterSite  = (clusterSize_ / 2) * bathSitesPerSite_;
		SizeType lastP1ClusterSite = firstClusterSite + clusterSize_;
		if (i >= firstClusterSite && i < lastP1ClusterSite)
			return PairType(-1, -1);

		SizeType middle = linSize_ / 2;
		SizeType cs     = clusterSize_ / 2;
		// now i is in the bath:
		if (i < middle) { // i is in the system
			return PairType(i % cs + firstClusterSite, i / cs);
		}
		// is in the bath and in the environ:
		SizeType iprime = i - lastP1ClusterSite;
		SizeType offset = lastP1ClusterSite - cs;
		return PairType(iprime % cs + offset, iprime / cs);
	}

	// assumes i1 and i2 are in the cluster
	// if connected return true, else false
	bool connectedInCluster(SizeType i1, SizeType i2) const
	{
		ladderize(i1, i2);
		return ladder_->connected(i1, i2);
	}

	// assumes i1 and i2 are connected and in the cluster
	SizeType calcDirInCluster(SizeType i1, SizeType i2) const
	{
		ladderize(i1, i2);
		return ladder_->calcDir(i1, i2);
	}

	// assumes i1 and i2 are in the cluster
	bool fringeInCluster(SizeType i, SizeType smax, SizeType emin) const
	{
		SizeType firstClusterSite = (clusterSize_ / 2) * bathSitesPerSite_;
		i -= firstClusterSite;
		smax -= firstClusterSite;
		emin -= firstClusterSite;
		return ladder_->fringe(i, smax, emin);
	}

	// assumes i1 and i2 are connected and in the cluster
	SizeType handleInCluster(SizeType i1, SizeType i2) const
	{
		ladderize(i1, i2);
		return ladder_->handle(i1, i2);
	}

	void ladderize(SizeType& i1, SizeType& i2) const
	{
		SizeType firstClusterSite = (clusterSize_ / 2) * bathSitesPerSite_;
		i1 -= firstClusterSite;
		i2 -= firstClusterSite;
	}

	SizeType    linSize_;
	SizeType    bathSitesPerSite_;
	SizeType    clusterSize_;
	LadderType* ladder_;
}; // class LadderBath
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_H
