// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file GeometryDirection.h
 *
 *  FIXME DOC. TBW
 */
#ifndef GEOMETRY_DIR_H
#define GEOMETRY_DIR_H
#include <PsimagLite/Io/IoSerializerStub.h>
#include <PsimagLite/Matrix.h>
#include <cassert>

namespace PsimagLite {

template <typename ComplexOrRealType, typename GeometryBaseType> class GeometryDirection {

	using MatrixType = Matrix<ComplexOrRealType>;
	using RealType   = typename Real<ComplexOrRealType>::Type;

public:

	enum InternalDofEnum
	{
		GENERAL,
		SPECIFIC
	};

	struct Auxiliary {

		Auxiliary(bool c, SizeType d, InternalDofEnum idof_, SizeType orbitals_)
		    : constantValues(c)
		    , dirId(d)
		    , idof(idof_)
		    , orbitals(orbitals_)
		{ }

		void write(PsimagLite::String label, IoSerializer& ioSerializer) const
		{
			ioSerializer.createGroup(label);
			ioSerializer.write(label + "/constantValues", constantValues);
			ioSerializer.write(label + "/dirId", dirId);
			ioSerializer.write(label + "/edof", idof);
		}

		bool            constantValues;
		SizeType        dirId;
		InternalDofEnum idof;
		SizeType        orbitals;
	}; // struct Auxiliary

	template <typename IoInputter>
	GeometryDirection(IoInputter&             io,
	                  const Auxiliary&        aux,
	                  const GeometryBaseType* geometryFactory)
	    : aux_(aux)
	    , geometryBase_(geometryFactory)
	{
		SizeType n = (aux.idof == GENERAL) ? 0 : getVectorSize();

		if (aux.idof == GENERAL) {
			geometryBase_->set(rawHoppings_, aux.orbitals);
			return;
		}

		assert(aux.idof == SPECIFIC);

		String connectors  = "Connectors";
		String savedPrefix = io.prefix();
		io.prefix() += "dir" + ttos(aux.dirId) + ":";

		if (aux.orbitals > 1) {
			if (n == 0)
				n = 1;
			for (SizeType i = 0; i < n; i++) {
				MatrixType m;
				String     extraString = (n > 1 && io.version() > 2) ? ttos(i) : "";
				io.read(m, connectors + extraString);
				dataMatrices_.push_back(m);
				if (aux.orbitals != m.rows() || aux.orbitals != m.cols())
					throw RuntimeError("Connectors must be matrices of "
					                   "rows=cols= "
					                   + ttos(aux.orbitals) + "\n");
			}
		} else {
			io.read(dataNumbers_, connectors);
			if (dataNumbers_.size() != n) {
				String s(__FILE__);
				s += " " + ttos(dataNumbers_.size()) + " != " + ttos(n) + "\n";
				throw RuntimeError(s.c_str());
			}
		}

		io.prefix() = savedPrefix;
	}

	void write(PsimagLite::String label, IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		aux_.write(label + "/aux_", ioSerializer);
		// geometryBase_
		ioSerializer.write(label + "/dataNumbers_", dataNumbers_);
		ioSerializer.write(label + "/dataMatrices_", dataMatrices_);
		rawHoppings_.write(label + "/rawHoppings_", ioSerializer);
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, String) const
	{
		return 0;
	}

	ComplexOrRealType operator()(SizeType i, SizeType edof1, SizeType j, SizeType edof2) const
	{
		SizeType ii = edof1 + i * aux_.orbitals;
		SizeType jj = edof2 + j * aux_.orbitals;
		if (aux_.idof == GENERAL) {
			assert(edof1 < aux_.orbitals && edof2 < aux_.orbitals);
			return rawHoppings_(ii, jj);
		}

		assert(aux_.idof == SPECIFIC);

		SizeType h = (constantValues()) ? 0 : geometryBase_->handle(i, j);

		bool isMatrix = (aux_.orbitals > 1);
		if (!isMatrix) {
			assert(dataNumbers_.size() > h);
			return dataNumbers_[h];
		}

		if (ii < jj)
			return PsimagLite::conj(operator()(j, edof2, i, edof1));

		assert(dataMatrices_.size() > h);

		bool b = (dataMatrices_[h].rows() > edof1 && dataMatrices_[h].cols() > edof2);

		assert(b || (dataMatrices_[h].rows() > edof2 && dataMatrices_[h].cols() > edof1));

		ComplexOrRealType tmp
		    = (b) ? dataMatrices_[h](edof1, edof2) : dataMatrices_[h](edof2, edof1);
		int signChange = geometryBase_->signChange(i, j);
		return tmp * static_cast<RealType>(signChange);
	}

	bool constantValues() const { return aux_.constantValues; }

	friend std::ostream& operator<<(std::ostream& os, const Auxiliary& a)
	{
		os << "constantValues=" << a.constantValues << "\n";
		os << "dirId=" << a.dirId << "\n";
		os << "edof=" << a.idof << "\n";

		return os;
	}

	friend std::ostream& operator<<(std::ostream& os, const GeometryDirection& gd)
	{
		os << "#GeometryDirectionAuxiliary\n";
		os << gd.aux_;

		bool isMatrix = (gd.aux_.orbitals > 1 || gd.aux_.idof == SPECIFIC);

		if (!isMatrix) {
			os << "#GeometryNumbersSize=" << gd.dataNumbers_.size() << "\n";
			os << "#GeometryNumbers=";
			for (SizeType i = 0; i < gd.dataNumbers_.size(); i++) {
				os << gd.dataNumbers_[i] << " ";
			}
			os << "\n";
		} else {
			os << "#GeometryMatrixSize=" << gd.dataMatrices_.size() << "\n";
			for (SizeType i = 0; i < gd.dataMatrices_.size(); i++)
				os << gd.dataMatrices_[i];
		}
		return os;
	}

private:

	SizeType getVectorSize()
	{
		return (aux_.constantValues) ? 1 : geometryBase_->getVectorSize(aux_.dirId);
	}

	Auxiliary                                aux_;
	const GeometryBaseType*                  geometryBase_;
	typename Vector<ComplexOrRealType>::Type dataNumbers_;
	typename Vector<MatrixType>::Type        dataMatrices_;
	MatrixType                               rawHoppings_;
}; // class GeometryDirection
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_DIR_H
