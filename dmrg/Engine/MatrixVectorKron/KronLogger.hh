#ifndef KRONLOGGER_HH
#define KRONLOGGER_HH
#include "Basis.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "MatrixMarket.hh"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <optional>

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
template <typename InitKronType_> class KronLogger {
public:

	using InitKronType            = InitKronType_;
	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using ComplexOrRealType       = typename InitKronType::ComplexOrRealType;
	using MatrixMarketType        = MatrixMarket<ComplexOrRealType>;

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
	void logVector(const std::vector<ComplexOrRealType>& v)
	{
		if (!fout_.has_value() || !active_)
			return;

		std::ios_base::fmtflags flags(fout_->flags());
		*fout_ << std::defaultfloat << std::setprecision(precision<ComplexOrRealType>());

		*fout_ << "Vector\n";
		*fout_ << v;

		fout_->flags(flags);
	}

	/*!
	 * \brief Logs the innermost ("connections") loop of the KronDiag algorithm
	 *
	 * \param[in] outPatch The outer patch index
	 * \param[in] inPatch The inner patch index
	 * \param[in] ic The connection index
	 */
	void logHamiltonian(SizeType outPatch, SizeType inPatch, SizeType ic, bool applyAdjoint)
	{
		if (!fout_.has_value() || !active_)
			return;

		std::ios_base::fmtflags flags(fout_->flags());
		*fout_ << std::defaultfloat << std::setprecision(precision<ComplexOrRealType>());

		const ArrayOfMatStructType& xiStruct = init_kron_.xc(ic);
		const ArrayOfMatStructType& yiStruct = init_kron_.yc(ic);

		const bool performTranspose
		    = (!applyAdjoint && init_kron_.useLowerPart() && (outPatch < inPatch));
		const SizeType sourceOutPatch = applyAdjoint ? inPatch : outPatch;
		const SizeType sourceInPatch  = applyAdjoint ? outPatch : inPatch;

		const MatrixDenseOrSparseType* Amat = performTranspose
		    ? xiStruct(inPatch, outPatch)
		    : xiStruct(sourceOutPatch, sourceInPatch);

		const MatrixDenseOrSparseType* Bmat = performTranspose
		    ? yiStruct(inPatch, outPatch)
		    : yiStruct(sourceOutPatch, sourceInPatch);

		*fout_
		    << "\n# nc, total, outPatch, inPatch, ic, applyAdjoint, performTranspose, offsetX, offsetY\n";
		*fout_ << init_kron_.connections() << ' '
		       << init_kron_.numberOfPatches(InitKronType::OLD) << ' ' << outPatch << ' '
		       << inPatch << ' ' << ic << ' ' << applyAdjoint << ' ' << performTranspose
		       << ' ' << init_kron_.offsetForPatches(InitKronType::NEW, outPatch) << ' '
		       << init_kron_.offsetForPatches(InitKronType::OLD, inPatch) << "\n";

		assert(Amat);
		*fout_ << "## A: " + matrixFormat(Amat->isDense()) + "\n";
		printMatrixDenseOrSparse(*Amat);

		assert(Bmat);
		*fout_ << "## B: " + matrixFormat(Bmat->isDense()) + "\n";
		printMatrixDenseOrSparse(*Bmat);

		fout_->flags(flags);
	}

	void sync() { active_ = false; }

private:

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
			tmp += ' ';
		}

		return tmp;
	}

	template <typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
	int precision() const
	{
		return std::numeric_limits<T>::max_digits10;
	}
	template <typename T, std::enable_if_t<!std::is_floating_point<T>::value, int> = 0>
	int precision() const
	{
		return std::numeric_limits<typename T::value_type>::max_digits10;
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
		size_t dot_index = filename.find_last_of('.');

		std::string root = filename.substr(0, dot_index);
		return "kron_" + root + "_" + ttos(n) + ".txt";
	}

	static SizeType               counter_;
	PsimagLite::ProgressIndicator progress_;
	std::optional<std::ofstream>  fout_;
	const InitKronType&           init_kron_;
	bool                          active_;
};

template <typename InitKronType> SizeType KronLogger<InitKronType>::counter_ = 0;
}
#endif // KRONLOGGER_HH
