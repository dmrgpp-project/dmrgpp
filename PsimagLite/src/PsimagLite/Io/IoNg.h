// SPDX-FileCopyrightText: Copyright (c) 2009-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file IoNg.h
 *
 *  This class handles Input/Output for PsimagLite
 */

#ifndef PSI_IO_NG_H
#define PSI_IO_NG_H

#include "PsimagLite/AllocatorCpu.h"
#include "PsimagLite/Io/IoNgSerializer.h"
#include "PsimagLite/Map.h"
#include "PsimagLite/Matrix.h"
#include "PsimagLite/Stack.h"
#include "PsimagLite/Vector.h"

#include <H5Cpp.h>
#include <iostream>
#include <typeinfo>

namespace PsimagLite {

/* PSIDOC IsRootUnDelegated
Root-undelegated types are one of the following.
\begin{lstlisting}
PSIDOCCOPY IsRootUnDelegatedCode
\end{lstlisting}
Root-undelegateds are always at least partially written to by \code{IoNg},
and, if needed, parts of it are delegated.
For example, all native types are written by \code{IoNg} directly into a single
dataset. \code{std::complex<T>} where \code{T} is a native type is written
directly by doubling the size of the array into a single dataset.
*/
/* PSIDOC_CODE_START IsRootUnDelegatedCode */
template <typename T> struct IsRootUnDelegated {
	enum
	{
		True = Loki::TypeTraits<T>::isArith || IsVectorLike<T>::True || IsStackLike<T>::True
		    || IsPairLike<T>::True || std::is_enum<T>::value || IsEnumClass<T>::value
		    || IsStringLike<T>::True
	};
};
/* PSIDOC_CODE_END */

class IoNg {

public:

	/*
	        H5F_ACC_TRUNC - Truncate file, if it already exists,
	        erasing all data previously stored in the file.
	        H5F_ACC_EXCL - Fail if file already exists. H5F_ACC_TRUNC
	        and H5F_ACC_EXCL are mutually exclusive
	        H5F_ACC_RDONLY - Open file as read-only, if it already exists,
	   and fail, otherwise H5F_ACC_RDWR - Open file for read/write, if it
	   already exists, and fail, otherwise
	*/
	enum OpenMode
	{
		ACC_TRUNC,
		ACC_EXCL,
		ACC_RDONLY,
		ACC_RDW
	};

	static void dontPrintDebug() { IoNgSerializer::dontPrintDebug(); }

	class Out {

	public:

		using Serializer       = IoNgSerializer;
		using VectorStringType = std::vector<String>;

		Out(const String& filename, OpenMode mode)
		    : ioNgSerializer_(filename, modeToH5(mode))
		{ }

		void flush() { ioNgSerializer_.flush(); }

		const String& filename() const { return ioNgSerializer_.filename(); }

		void open(String filename, OpenMode mode)
		{
			ioNgSerializer_.open(filename, modeToH5(mode));
		}

		void close() { ioNgSerializer_.close(); }

		void createGroup(String groupName) { ioNgSerializer_.createGroup(groupName); }

		template <typename T>
		void writeVectorEntry(T x, PsimagLite::String str, SizeType counter)
		{
			if (counter == 0)
				createGroup(str);

			ioNgSerializer_.write(str + "/" + ttos(counter), x);
			ioNgSerializer_.write(str + "/Size",
			                      counter + 1,
			                      (counter == 0) ? IoNgSerializer::NO_OVERWRITE
			                                     : IoNgSerializer::ALLOW_OVERWRITE);
		}

		template <typename T>
		void write(std::stack<T>&            what,
		           String                    name2,
		           IoNgSerializer::WriteMode mode = IoNgSerializer::NO_OVERWRITE,
		           typename EnableIf<!IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			ioNgSerializer_.write(name2, what, mode);
		}

		template <typename T>
		void write(const T&                  what,
		           String                    name2,
		           IoNgSerializer::WriteMode mode = IoNgSerializer::NO_OVERWRITE,
		           typename EnableIf<IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			ioNgSerializer_.write(name2, what, mode);
		}

		template <typename T>
		void write(const T& what,
		           String   name2,
		           typename EnableIf<!IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			what.write(name2, ioNgSerializer_);
		}

		template <typename T>
		void write(const T&                  what,
		           String                    name2,
		           IoNgSerializer::WriteMode mode,
		           typename EnableIf<!IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			what.write(name2, ioNgSerializer_, mode);
		}

		template <typename T>
		void overwrite(const T& what,
		               String   name2,
		               typename EnableIf<IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			ioNgSerializer_.overwrite(name2, what);
		}

		template <typename T>
		void overwrite(const T& what,
		               String   name2,
		               typename EnableIf<!IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			what.overwrite(name2, ioNgSerializer_);
		}

		IoNgSerializer& serializer() { return ioNgSerializer_; }

	private:

		Out(const Out&);

		Out& operator=(const Out&);

		static unsigned int modeToH5(OpenMode mode)
		{
			switch (mode) {
			case ACC_TRUNC:
				return H5F_ACC_TRUNC;
			case ACC_EXCL:
				return H5F_ACC_EXCL;
			case ACC_RDONLY:
				return H5F_ACC_RDONLY;
			case ACC_RDW:
				return H5F_ACC_RDWR;
			}

			throw RuntimeError("IoNg:: wrong open mode\n");
		}

		IoNgSerializer ioNgSerializer_;
	};

	class In {

	public:

		using LongIntegerType = int long;
		using LongSizeType    = unsigned int long;

		In(String filename)
		    : ioNgSerializer_(filename, H5F_ACC_RDONLY)
		{ }

		const String& filename() const { return ioNgSerializer_.filename(); }

		void open(String filename) { ioNgSerializer_.open(filename, H5F_ACC_RDONLY); }

		void close() { ioNgSerializer_.close(); }

		template <typename SomeType> void readLastVectorEntry(SomeType& x, String s)
		{
			int total = 0;
			ioNgSerializer_.read(total, s + "/Size");

			if (total <= 0)
				throw RuntimeError("Error reading last instance of " + s + "\n");

			ioNgSerializer_.read(x, s + "/" + ttos(--total));
		}

		template <typename T>
		void read(T&     what,
		          String name,
		          typename EnableIf<IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			ioNgSerializer_.read(what, name);
		}

		template <typename T>
		void read(T&     what,
		          String name,
		          typename EnableIf<!IsRootUnDelegated<T>::True, int*>::Type = 0)
		{
			what.read(name, ioNgSerializer_);
		}

		IoNgSerializer& serializer() { return ioNgSerializer_; }

	private:

		In(const In&);

		In& operator=(const In&);

		IoNgSerializer ioNgSerializer_;
	};
}; // class IoNg

template <> struct IsInputLike<IoNg::In> {
	enum
	{
		True = true
	};
};

template <> struct IsOutputLike<IoNg::Out> {
	enum
	{
		True = true
	};
};

} // namespace PsimagLite

/*@}*/
#endif
