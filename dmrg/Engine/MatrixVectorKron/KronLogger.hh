#ifndef KRONLOGGER_HH
#define KRONLOGGER_HH
#include "ApplyOperatorLocal.h"
#include "Basis.h"
#include "BasisWithOperators.h"
#include "FermionSign.h"
#include "InitKronHamiltonian.h"
#include "LeftRightSuper.h"
#include "MatrixMarket.hh"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <limits>
#include <optional>
#include <vector>

namespace Dmrg {

/*!
 * The Kron logger logs data that is sent to KronUtil/MatrixDenseOrSparse.h
 * function kronMult()
 *
 * It logs first the constant vector v in H*v, then it logs data in three
 * stages, as you can see in member functions one(), two(), and three().
 *
 * DMRG++ calls H*v or matrixVector(x, y) many times so the user
 * enters a start via KroneckerDumperBegin, and an end via KroneckerDumperEnd,
 * and also enables KroneckerDumper in SolverOptions; see file input31.ain
 * in the TestSuite.
 *
 * If the count for H*v is greater or equal start and less than end, then
 * the logging is active, and each dataset going to kronMult is logged in
 * a separate file, with a counter at the end and with a rootname
 * based on the output file of the run; see buildFilename() below.
 */
template <typename ModelType> class KronLogger {
public:

	using InitKronType            = InitKronHamiltonian<ModelType>;
	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using ComplexOrRealType       = typename InitKronType::ComplexOrRealType;
	using SparseMatrixType        = typename InitKronType::SparseMatrixType;
	using LeftRightSuperType      = typename InitKronType::LeftRightSuperType;
	using VectorType              = typename InitKronType::VectorType;
	using RealType                = typename ModelType::RealType;
	using VectorWithOffsetType    = VectorWithOffsets<ComplexOrRealType>;
	using ApplyOperatorType = ApplyOperatorLocal<LeftRightSuperType, VectorWithOffsetType>;
	using MatrixMarketType  = MatrixMarket<ComplexOrRealType>;

	/*!
	 * \brief CONSTRUCTOR
	 *
	 * \param[in] init_kron The InitKronHamiltonian object as const reference
	 *
	 */
	KronLogger(const InitKronType& init_kron)
	    : progress_("KronLogger")
	    , fout_()
	    , init_kron_(init_kron)
	    , active_(true)
	{
		if (!init_kron.params().options.isSet("KroneckerDumper")) {
			return;
		}

		const std::string message
		    = "KronLogger: Hello from ctor, counter_=" + ttos(counter_);
		progress_.printline(message, std::cout);

		++counter_;

		// Here counter_ has a value larger-by-one to what was printed
		// so we use counter_ - 1 instead to compare with start and end
		// it includes start but not end
		// assumes start == end == 0 or start < end which where
		// checked in the ParametersForDmrg reading

		SizeType end = init_kron.params().dumperEnd;

		if (end == 0 || counter_ >= end + 1) {
			return;
		}

		SizeType start = init_kron.params().dumperBegin;
		assert((start == 0 && end == 0) || (start < end));

		if (counter_ < start + 1) {
			return;
		}

		assert(counter_ >= start + 1);
		std::string filename
		    = buildFilename(init_kron.params().filename, counter_ - start - 1);
		fout_.emplace(filename);
		if (!fout_ or !fout_->good() or fout_->bad()) {
			fout_.reset();
			err(std::string("Failed to create KronLogger file ") + filename + "\n");
		}

		printMetadata(message);
		buildLocalOperator();
	}

	/*!
	 * \brief DESTRUCTOR
	 */
	~KronLogger()
	{
		if (!fout_.has_value())
			return;

		printMetadata("KronLogger: Bye from dtor");
	}

	/*!
	 * \brief Logs the constant vector v for H*v
	 *
	 * \param[in] v The vector
	 */
	void vector(const std::vector<ComplexOrRealType>& v)
	{
		if (!fout_.has_value() || !active_)
			return;

		*fout_ << "Vector\n";
		const std::streamsize oldPrecision = fout_->precision();
		if (!init_kron_.params().dumperOperator.empty())
			fout_->precision(std::numeric_limits<RealType>::max_digits10);
		*fout_ << v;
		if (!init_kron_.params().dumperOperator.empty()) {
			const ComplexOrRealType expectation = verifyLocalOperator(v);
			*fout_ << "OperatorExpectation=" << PsimagLite::real(expectation) << " "
			       << PsimagLite::imag(expectation) << "\n";
		}
		fout_->precision(oldPrecision);
	}

	/*!
	 * \brief Logs the outer loop of the KronDiag algorithm
	 *
	 * \param[in] outPatch The outer patch index
	 */
	void one(SizeType outPatch)
	{
		if (!fout_.has_value() || !active_)
			return;

		SizeType nC      = init_kron_.connections();
		SizeType total   = init_kron_.numberOfPatches(InitKronType::OLD);
		SizeType offsetX = init_kron_.offsetForPatches(InitKronType::NEW, outPatch);
		*fout_ << "outPatch=" << outPatch << "\n";
		*fout_ << "nC=" << nC << " total=" << total << " offsetX=" << offsetX << "\n";
	}

	/*!
	 * \brief Logs the middle loop of the KronDiag algorithm
	 *
	 * \param[in] inPatch The inner patch index
	 */
	void two(SizeType inPatch)
	{
		if (!fout_.has_value() || !active_)
			return;

		SizeType offsetY = init_kron_.offsetForPatches(InitKronType::OLD, inPatch);
		*fout_ << separationLevel(1) << "inPatch=" << inPatch << "\n";
		*fout_ << separationLevel(1) << "offsetY=" << offsetY << "\n";
	}

	/*!
	 * \brief Logs the innermost ("connections") loop of the KronDiag algorithm
	 *
	 * \param[in] outPatch The outer patch index
	 * \param[in] inPatch The inner patch index
	 * \param[in] ic The connection index
	 */
	void three(SizeType outPatch, SizeType inPatch, SizeType ic)
	{
		if (!fout_.has_value() || !active_)
			return;

		const ArrayOfMatStructType& xiStruct = init_kron_.xc(ic);
		const ArrayOfMatStructType& yiStruct = init_kron_.yc(ic);

		const bool performTranspose = (init_kron_.useLowerPart() && (outPatch < inPatch));

		const MatrixDenseOrSparseType* Amat
		    = performTranspose ? xiStruct(inPatch, outPatch) : xiStruct(outPatch, inPatch);

		const MatrixDenseOrSparseType* Bmat
		    = performTranspose ? yiStruct(inPatch, outPatch) : yiStruct(outPatch, inPatch);

		*fout_ << separationLevel(2) << "ic=" << ic << "\n";
		*fout_ << separationLevel(2) << "performTranspose=" << performTranspose << "\n";

		assert(Amat);
		*fout_ << separationLevel(2)
		       << "Matrix A follows in format: " + matrixFormat(Amat->isDense()) + "\n";
		printMatrixDenseOrSparse(*Amat);

		assert(Bmat);
		*fout_ << separationLevel(2)
		       << "Matrix B follows in format: " + matrixFormat(Bmat->isDense()) + "\n";
		printMatrixDenseOrSparse(*Bmat);
	}

	void sync() { active_ = false; }

private:

	struct MatrixEntry {

		SizeType          row;
		SizeType          col;
		ComplexOrRealType value;
	};

	void buildLocalOperator()
	{
		if (!fout_.has_value() || init_kron_.params().dumperOperator.empty())
			return;

		const LeftRightSuperType& lrs       = init_kron_.lrs(InitKronType::NEW);
		const auto                direction = init_kron_.direction();
		assert(!lrs.right().block().empty());
		localOperatorSite_ = (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		    ? lrs.right().block()[0] - 1
		    : lrs.right().block()[0];

		const ModelType& model         = init_kron_.model();
		const SizeType   numberOfSites = model.superGeometry().numberOfSites();
		if (localOperatorSite_ >= numberOfSites)
			err("KronLogger: inferred local-operator site is outside the geometry\n");
		const auto& localOperator = model.naturalOperator(
		    init_kron_.params().dumperOperator, localOperatorSite_, 0);
		typename PsimagLite::Vector<bool>::Type oddElectrons;
		model.findOddElectronsOfOneSite(oddElectrons, localOperatorSite_);
		FermionSign       fermionSign(lrs.left(), oddElectrons);
		ApplyOperatorType applyOperator(lrs, false);

		const SizeType sectorSize   = init_kron_.size(InitKronType::NEW);
		const SizeType sectorOffset = init_kron_.offset(InitKronType::NEW);
		const SizeType sector       = lrs.super().findPartitionNumber(sectorOffset);
		const SizeType splitSize    = model.hilbertSize(localOperatorSite_);
		if (localOperator.getCRS().rows() != splitSize
		    || localOperator.getCRS().cols() != splitSize)
			err("KronLogger: local-operator dimensions do not match the site Hilbert space\n");
		VectorType               kronBasis(sectorSize, 0);
		VectorType               sectorBasis;
		VectorType               sectorResult(sectorSize, 0);
		VectorType               kronResult;
		std::vector<MatrixEntry> entries;
		entries.reserve(sectorSize);
		const auto border
		    = (localOperatorSite_ == 0 || localOperatorSite_ + 1 == numberOfSites)
		    ? ApplyOperatorType::BORDER_YES
		    : ApplyOperatorType::BORDER_NO;

		for (SizeType col = 0; col < sectorSize; ++col) {
			kronBasis[col] = 1;
			init_kron_.kronToSector(sectorBasis, kronBasis);

			VectorWithOffsetType source;
			source.set(sectorBasis, sector, lrs.super());
			VectorWithOffsetType destination;
			applyOperator(destination,
			              source,
			              localOperator,
			              fermionSign,
			              splitSize,
			              direction,
			              border);

			destination.extract(sectorResult, sector);
			if (sectorResult.empty())
				sectorResult.resize(sectorSize, 0);
			assert(sectorResult.size() == sectorSize);
			init_kron_.sectorToKron(kronResult, sectorResult);
			for (SizeType row = 0; row < sectorSize; ++row) {
				if (kronResult[row] != ComplexOrRealType(0))
					entries.push_back({ row, col, kronResult[row] });
			}

			kronBasis[col] = 0;
			for (SizeType i = 0; i < sectorSize; ++i)
				sectorResult[i] = 0;
		}

		std::sort(entries.begin(),
		          entries.end(),
		          [](const MatrixEntry& a, const MatrixEntry& b)
		          { return (a.row < b.row || (a.row == b.row && a.col < b.col)); });
		localOperator_.resize(sectorSize, sectorSize);
		SizeType entry = 0;
		for (SizeType row = 0; row < sectorSize; ++row) {
			localOperator_.setRow(row, entry);
			while (entry < entries.size() && entries[entry].row == row) {
				localOperator_.pushCol(entries[entry].col);
				localOperator_.pushValue(entries[entry].value);
				++entry;
			}
		}

		localOperator_.setRow(sectorSize, entry);
		localOperator_.checkValidity();
		printLocalOperator(direction, border);
		progress_.printline("KronLogger: built " + init_kron_.params().dumperOperator
		                        + " at site " + ttos(localOperatorSite_) + " as "
		                        + ttos(sectorSize) + "x" + ttos(sectorSize) + " with "
		                        + ttos(localOperator_.nonZeros()) + " nonzeros",
		                    std::cout);
	}

	ComplexOrRealType verifyLocalOperator(const VectorType& kronVector)
	{
		const SizeType sectorSize = init_kron_.size(InitKronType::NEW);
		if (sectorSize == 0 || kronVector.size() != sectorSize
		    || localOperator_.rows() != sectorSize || localOperator_.cols() != sectorSize)
			err("KronLogger: local-operator verification dimension mismatch\n");

		const LeftRightSuperType& lrs           = init_kron_.lrs(InitKronType::NEW);
		const auto                direction     = init_kron_.direction();
		const ModelType&          model         = init_kron_.model();
		const auto&               localOperator = model.naturalOperator(
                    init_kron_.params().dumperOperator, localOperatorSite_, 0);
		typename PsimagLite::Vector<bool>::Type oddElectrons;
		model.findOddElectronsOfOneSite(oddElectrons, localOperatorSite_);
		FermionSign       fermionSign(lrs.left(), oddElectrons);
		ApplyOperatorType applyOperator(lrs, false);
		const SizeType    sectorOffset  = init_kron_.offset(InitKronType::NEW);
		const SizeType    sector        = lrs.super().findPartitionNumber(sectorOffset);
		const SizeType    splitSize     = model.hilbertSize(localOperatorSite_);
		const SizeType    numberOfSites = model.superGeometry().numberOfSites();
		const auto        border
		    = (localOperatorSite_ == 0 || localOperatorSite_ + 1 == numberOfSites)
		    ? ApplyOperatorType::BORDER_YES
		    : ApplyOperatorType::BORDER_NO;

		const auto applyDirect = [&](const VectorType& input)
		{
			VectorType sectorVector;
			init_kron_.kronToSector(sectorVector, input);
			VectorWithOffsetType source;
			source.set(sectorVector, sector, lrs.super());
			VectorWithOffsetType destination;
			applyOperator(destination,
			              source,
			              localOperator,
			              fermionSign,
			              splitSize,
			              direction,
			              border);

			VectorType directSectorResult(sectorSize, 0);
			destination.extract(directSectorResult, sector);
			if (directSectorResult.empty())
				directSectorResult.resize(sectorSize, 0);
			VectorType directKronResult;
			init_kron_.sectorToKron(directKronResult, directSectorResult);
			return directKronResult;
		};

		const VectorType      directKronResult = applyDirect(kronVector);
		std::vector<SizeType> sampleColumns    = { 0, sectorSize / 2, sectorSize - 1 };
		std::sort(sampleColumns.begin(), sampleColumns.end());
		sampleColumns.erase(std::unique(sampleColumns.begin(), sampleColumns.end()),
		                    sampleColumns.end());
		for (const SizeType col : sampleColumns) {
			VectorType basisVector(sectorSize, 0);
			basisVector[col]              = 1;
			const VectorType directColumn = applyDirect(basisVector);
			*fout_ << "OperatorDirectColumn=" << col << "\n";
			*fout_ << directColumn;
		}

		return scalarProduct(kronVector, directKronResult);
	}

	void printLocalOperator(ProgramGlobals::DirectionEnum          direction,
	                        typename ApplyOperatorType::BorderEnum border)
	{
		assert(fout_.has_value());
		*fout_ << "LocalOperator\n";
		*fout_ << "OperatorName=" << init_kron_.params().dumperOperator << "\n";
		*fout_ << "OperatorSite=" << localOperatorSite_ << "\n";
		*fout_ << "Rows=" << localOperator_.rows() << " Cols=" << localOperator_.cols()
		       << "\n";
		*fout_ << "BasisOrder=KronPatchPacked\n";
		*fout_ << "Direction=" << ProgramGlobals::toString(direction) << "\n";
		*fout_ << "Border=" << (border == ApplyOperatorType::BORDER_YES) << "\n";
		*fout_ << "MatrixAction=dest=O*src\n";
		*fout_ << "MatrixFormat=matrix-market-coordinate\n";
		const std::streamsize oldPrecision = fout_->precision();
		fout_->precision(std::numeric_limits<RealType>::max_digits10);
		MatrixMarketType(localOperator_).print(*fout_);
		fout_->precision(oldPrecision);
	}

	void printMetadata(const std::string& message)
	{
		if (PsimagLite::Concurrency::rank() != 0)
			return;

		*fout_ << "## ";
		progress_.printline(message, *fout_);
	}

	static std::string separationLevel(SizeType n)
	{
		std::string tmp;
		for (SizeType i = 0; i < n; ++i) {
			tmp += " ";
		}

		return tmp;
	}

	static std::string matrixFormat(bool b) { return (b) ? "dense" : "sparse"; }

	void printMatrixDenseOrSparse(const MatrixDenseOrSparseType& mat)
	{
		if (!fout_) {
			err("printMatrixDenseOrSparse: InternalError: KronLogger not enabled!?\n");
		}

		if (mat.isDense()) {
			*fout_ << mat.dense();
		} else {
			MatrixMarketType matrix_market(mat.sparse());
			matrix_market.print(*fout_);
		}
	}

	static std::string buildFilename(const std::string& filename, SizeType n)
	{
		size_t dot_index = filename.find_last_of(".");

		std::string root = filename.substr(0, dot_index);
		return "kron_" + root + "_" + ttos(n) + ".txt";
	}

	static SizeType               counter_;
	PsimagLite::ProgressIndicator progress_;
	std::optional<std::ofstream>  fout_;
	const InitKronType&           init_kron_;
	bool                          active_;
	SparseMatrixType              localOperator_;
	SizeType                      localOperatorSite_ = 0;
};

template <typename ModelType> SizeType KronLogger<ModelType>::counter_ = 0;
}
#endif // KRONLOGGER_HH
