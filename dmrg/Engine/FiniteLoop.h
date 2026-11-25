#ifndef DMRG_FINITELOOP_H
#define DMRG_FINITELOOP_H

#include "ProgramGlobals.h"
#include "TruncationControl.h"
#include <PsimagLite/Io/IoSerializerStub.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>
#include <array>
#include <iostream>
#include <utility>

namespace Dmrg {
/* PSIDOC FiniteLoop
\subsection{Finite loops}

To enable finite loops, make sure that `nofiniteloops` is \emph{not} present
under `SolverOptions=`. In Ainur input, `FiniteLoops` is a matrix with one row
per finite movement:

\begin{verbatim}
FiniteLoops = [[7, 200, 0], [-7, 200, @save]];
\end{verbatim}

For each row, the first entry is the number of sites to move forward if
positive or backward if negative, and the second entry is the number of kept
states ($m$). The third entry controls per-sweep options. The first movement
starts where the infinite loop ends, at the middle of the lattice.

\subsection{The third entry in the triplet}
See \texttt{doc/FiniteLoopThirdEntry.md} for the current syntax, options,
constraints, and examples.
\subsection{Caveats and Troubleshooting}

If `nofiniteloops` is an option in the options line of the input file, then
the \verb!FiniteLoops! line is ignored and no finite loops are done. In this
case, DMRG++ stops when the infinite algorithm has finished.

Make sure you don't fall off the lattice, by going forward or backwards too much.
Remember that at least one site must remain for the system part of the lattice.
So on a 16 site chain, when you start the finite loops you're at the middle, you
can go forward at most 7 sites, and backwards at most 7 sites.

There is some checking done to the finite loops input,
but you might find that it's not comprehensive.
*/
template <typename RealType> class FiniteLoop {

public:

	using VectorStringType      = PsimagLite::Vector<PsimagLite::String>::Type;
	using TruncationControlType = TruncationControl<RealType>;

	FiniteLoop(int                          sl,
	           SizeType                     ks,
	           PsimagLite::String           str,
	           const TruncationControlType& truncationControl)
	    : stepLength_(sl)
	    , keptStates_(ks)
	    , bitField_(0)
	    , truncationControl_(truncationControl)
	{
		str = removeCharsAtStart(str, "\"");
		str = removeCharsAtEnd(str, "\"");

		VectorStringType tokens;
		PsimagLite::split(tokens, str, ":");
		bool hasBitField = false;
		bool hasAtField  = false;
		for (SizeType i = 0; i < tokens.size(); ++i)
			procToken(hasBitField, hasAtField, tokens[i]);

		warnIfFiniteMlessThanMin();
		checkBitField();
	}

	void write(PsimagLite::String                    label,
	           PsimagLite::IoSerializer&             ioSerializer,
	           PsimagLite::IoNgSerializer::WriteMode mode
	           = PsimagLite::IoNgSerializer::NO_OVERWRITE) const
	{
		PsimagLite::String root = label;

		if (mode != PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
			ioSerializer.createGroup(root);
		ioSerializer.write(root + "/stepLength", stepLength_, mode);
		ioSerializer.write(root + "/keptStates", keptStates_, mode);
		ioSerializer.write(root + "/saveOption", bitField_, mode);
	}

	int stepLength() const { return stepLength_; }

	SizeType keptStates() const { return keptStates_; }

	const TruncationControlType& truncationControl() const { return truncationControl_; }

	bool wants(PsimagLite::String what) const
	{
		PsimagLite::String whatLower = ProgramGlobals::toLower(what);
		const OptionType*  option    = findOption(whatLower);
		assert(option != nullptr);
		return (bitField_ & option->second);
	}

	void print(std::ostream& os) const
	{
		os << stepLength_ << " " << keptStates_ << " " << bitField_;
	}

private:

	using OptionType  = std::pair<const char*, SizeType>;
	using OptionsType = std::array<OptionType, 6>;

	static const OptionsType& options()
	{
		static const OptionsType values = { { { "save", 1 },
			                              { "onlyfastwft", 2 },
			                              { "onlyslowwft", 4 },
			                              { "randomguess", 8 },
			                              { "multisitepush", 16 },
			                              { "onesitetruncation", 32 } } };
		return values;
	}

	static const OptionType* findOption(const PsimagLite::String& what)
	{
		const OptionsType& values = options();
		for (const OptionType& option : values) {
			if (what == option.first)
				return &option;
		}

		return nullptr;
	}

	// either @something or
	// setting=value
	void procToken(bool& hasBitField, bool& hasAtField, PsimagLite::String str)
	{
		if (isAdigit(str)) {
			dieIfTrue(hasBitField, "Bit field already set for this finite loop\n");
			dieIfTrue(hasAtField, "Bit field already set for this finite loop\n");
			hasBitField = true;
			bitField_   = PsimagLite::atoi(str);
			return;
		}

		if (str.length() < 2)
			err("FiniteLoop with " + str + " : syntax error (1)\n");

		if (str[0] == '@') {
			dieIfTrue(hasBitField, "Bit field already set for this finite loop\n");
			procAtValue(str.substr(1, str.length() - 1));
			hasAtField = true;
			return;
		}

		VectorStringType tokens;
		PsimagLite::split(tokens, str, "=");
		if (tokens.size() != 2)
			err("FiniteLoop with " + str + " : syntax error (2)\n");

		truncationControl_.setNameValue(tokens[0], tokens[1]);
	}

	void procAtValue(PsimagLite::String what)
	{
		const OptionType* option = findOption(what);
		if (option == nullptr) {
			err("Unknown finite-loop option @" + what
			    + ". Valid lowercase options are: @save, @onlyfastwft, @onlyslowwft, "
			      "@randomguess, @multisitepush, @onesitetruncation\n");
			return;
		}

		if (bitField_ & option->second)
			err("Finite-loop option @" + what + " already set\n");

		bitField_ |= option->second;
	}

	void warnIfFiniteMlessThanMin() const
	{
		if (keptStates_ >= truncationControl().mMin())
			return;

		std::cout << "WARNING: Triplet has m= " << keptStates_;
		std::cout << " which is less than its minimum m = " << truncationControl().mMin();
		std::cout << " as found in TruncationTolerance or set only for this loop\n";
	}

	void checkBitField() const
	{
		// Only 1 bit of bit 1, bit 2, and bit 3 may be set
		SizeType saveOption = bitField_;
		bool     flag       = false;
		for (SizeType i = 0; i < 3; ++i) {
			saveOption >>= 1; // ATTENTION: Discards bit 0 when i === 0
			if (!(saveOption & 1))
				continue;
			if (flag) {
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg << "Triplet 3rd: Only one bit of bit 1, bit 2, and bit 3 may "
				       "be set";
				msg << " VALUE= " << bitField_ << "\n";
				err(msg.str());
			} else {
				flag = true;
			}
		}

		bool wantsOneSiteTrunc = (bitField_ & 32);
		if (wantsOneSiteTrunc && !(bitField_ & 8))
			err(PsimagLite::String(
			        "A finite loop that wants one site truncation must also want")
			    + " random guess.\nIn other words, if bit 5 is set so must bit 3.\n");
	}

	static void dieIfTrue(bool b, PsimagLite::String msg)
	{
		if (!b)
			return;
		err(msg);
	}

	static bool isAdigit(PsimagLite::String str)
	{
		if (str.empty())
			return false;

		for (SizeType i = 0; i < str.length(); ++i)
			if (!std::isdigit(str[i]))
				return false;
		return true;
	}

	static PsimagLite::String removeCharsAtStart(PsimagLite::String str,
	                                             PsimagLite::String chars)
	{
		size_t start = str.find_first_not_of(chars);
		return (start == PsimagLite::String::npos) ? "" : str.substr(start);
	}

	static PsimagLite::String removeCharsAtEnd(PsimagLite::String str, PsimagLite::String chars)
	{
		size_t end = str.find_last_not_of(chars);
		return (end == PsimagLite::String::npos) ? "" : str.substr(0, end + 1);
	}

	const int             stepLength_; // how much to go right (+) or left (-)
	const SizeType        keptStates_; // kept states
	SizeType              bitField_;
	TruncationControlType truncationControl_;
};

} // namespace Dmrg

#endif // DMRG_FINITELOOP_H
