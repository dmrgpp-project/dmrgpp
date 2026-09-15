// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file ProgressIndicator.h
 *
 *  This class handles output to a progress indicator (usually the terminal)
 */

#ifndef PROGRESS_INDICATOR_H
#define PROGRESS_INDICATOR_H

#include "Concurrency.h"
#include "MemoryUsage.h"
#include "TypeToString.h"
#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

namespace PsimagLite {

class ProgressIndicator {

	static MemoryUsage   musage_;
	static OstringStream buffer_;
	static bool          bufferActive_;

public:

	ProgressIndicator(String caller, SizeType threadId = 0)
	    : threadId_(threadId)
	    , rank_(0)
	{
		if (threadId_ != 0)
			return;

		caller_ = caller;
		rank_   = Concurrency::rank();
	}

	static void updateBuffer(int signal)
	{
		if (bufferActive_) {
			pid_t  p = getpid();
			String outName("buffer");
			outName += ttos(p);
			outName += ".txt";
			std::ofstream fout(outName.c_str());
			fout << buffer_().str() << "\n";
			fout.close();
			buffer_().str("");
		}

		bufferActive_ = !bufferActive_;

		String bufferActive = (bufferActive_) ? "active" : "inactive";
		std::cerr << "ProgressIndicator: signal " << signal << " received.";
		std::cerr << " buffer is now " << bufferActive << "\n";
	}

	template <typename SomeOutputType> void printline(const String& s, SomeOutputType& os) const
	{
		if (threadId_ != 0)
			return;
		if (rank_ != 0)
			return;
		prefix(os);
		os << s << "\n";

		if (!bufferActive_)
			return;

		prefix(buffer_);
		buffer_() << s << "\n";
	}

	void printline(OstringStream& s, std::ostream& os) const
	{
		if (threadId_ != 0)
			return;
		if (rank_ != 0)
			return;
		prefix(os);
		os << s().str() << "\n";
		s().seekp(std::ios_base::beg);

		if (!bufferActive_)
			return;

		prefix(buffer_);
		buffer_() << s().str() << "\n";
		s().seekp(std::ios_base::beg);
	}

	void print(const String& something, std::ostream& os) const
	{
		if (threadId_ != 0)
			return;
		if (rank_ != 0)
			return;
		prefix(os);
		os << something;

		if (!bufferActive_)
			return;

		prefix(buffer_);
		buffer_() << something;
	}

	void printMemoryUsage()
	{
		musage_.update();
		String        vmPeak = musage_.findEntry("VmPeak:");
		String        vmSize = musage_.findEntry("VmSize:");
		OstringStream msg(std::cout.precision());
		msg() << "Current virtual memory is " << vmSize << " maximum was " << vmPeak;
		printline(msg, std::cout);

		if (!bufferActive_)
			return;

		buffer_() << "Current virtual memory is " << vmSize << " maximum was " << vmPeak;
		printline(buffer_, std::cout);
	}

	static MemoryUsage::TimeHandle time() { return musage_.time(); }

private:

	template <typename SomeOutputStreamType> void prefix(SomeOutputStreamType& os) const
	{
		const MemoryUsage::TimeHandle t       = musage_.time();
		const double                  seconds = t.millis();
		const SizeType                prec    = os.precision(3);
		prefixHelper(os) << caller_ << " "
		                 << "[" << std::fixed << seconds << "]: ";
		os.precision(prec);
	}

	OstringStream::OstringStreamType& prefixHelper(OstringStream& os) const { return os(); }

	template <typename SomeOutputStreamType>
	SomeOutputStreamType& prefixHelper(SomeOutputStreamType& os) const
	{
		return os;
	}

	String   caller_;
	SizeType threadId_;
	SizeType rank_;
}; // ProgressIndicator

} // namespace PsimagLite

/*@}*/
#endif
