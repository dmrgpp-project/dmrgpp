// SPDX-FileCopyrightText: Copyright (c) 2009-2013-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file MemoryUsage.h
 *
 */

#ifndef MEMORY_USAGE_H
#define MEMORY_USAGE_H
#include <fstream>
#include <iostream>
#include <sys/time.h>

namespace PsimagLite {

class MemoryUsage {

	static const SizeType MY_MAX_LINE = 40240;

public:

	class TimeHandle {

	public:

		TimeHandle()
		{
			startTime_.tv_sec  = 0;
			startTime_.tv_usec = 0;
			gettimeofday(&startTime_, 0);
		}

		TimeHandle(time_t s, suseconds_t u)
		{
			startTime_.tv_sec  = s;
			startTime_.tv_usec = u;
		}

		TimeHandle operator-(const TimeHandle& other) const
		{
			time_t      s = startTime_.tv_sec - other.startTime_.tv_sec;
			suseconds_t u = startTime_.tv_usec - other.startTime_.tv_usec;
			return TimeHandle(s, u);
		}

		TimeHandle& operator+=(const TimeHandle& other)
		{
			startTime_.tv_sec += other.startTime_.tv_sec;
			startTime_.tv_usec += other.startTime_.tv_usec;
			return *this;
		}

		time_t seconds() const { return startTime_.tv_sec; }

		suseconds_t useconds() const { return startTime_.tv_usec; }

		double millis() const
		{
			const double tmp = startTime_.tv_usec / 1000. + startTime_.tv_sec * 1000.;
			return tmp / 1000.;
		}

	private:

		timeval startTime_;
	};

	MemoryUsage(const String& myself = "")
	    : data_("")
	    , myself_(myself)
	{
		if (myself_ == "")
			myself_ = "/proc/self/status";
		update();
	}

	void update()
	{
		std::ifstream ifp(myself_.c_str());
		if (!ifp || !ifp.good() || ifp.bad())
			return;
		char tmp[MY_MAX_LINE];
		data_ = "";
		while (!ifp.eof()) {
			ifp.getline(tmp, MY_MAX_LINE);
			data_ += String(tmp);
			data_ += String("\n");
		}

		ifp.close();
	}

	String findEntry(const String& label)
	{
		long unsigned int x = data_.find(label);
		if (x == String::npos)
			return "NOT_FOUND";

		x += label.length();
		long unsigned int y   = data_.find("\n", x);
		SizeType          len = y - x;
		if (y == String::npos)
			len = data_.length() - x;
		String s2        = data_.substr(x, len);
		x                = 0;
		const SizeType n = s2.length();
		for (SizeType i = 0; i < n; ++i) {
			++x;
			if (s2.at(i) == ' ' || s2.at(i) == '\t')
				continue;
			else
				break;
		}

		if (x > 0)
			--x;
		len = s2.length() - x;
		return s2.substr(x, len);
	}

	TimeHandle time() const
	{
		TimeHandle nowtime;
		return nowtime - startTime_;
	}

private:

	String     data_;
	String     myself_;
	TimeHandle startTime_;
}; // class MemoryUsage

} // namespace PsimagLite

/*@}*/
#endif // MEMORY_USAGE_H
