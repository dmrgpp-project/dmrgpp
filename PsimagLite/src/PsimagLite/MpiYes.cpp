// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/
// This file is meant to be included by Mpi.cpp ONLY
/*! \file MpiYes.cpp
 *
 */

#include "MpiYes.h"

namespace PsimagLite {

namespace MPI {

	CommType COMM_WORLD = MPI_COMM_WORLD;

	template <> const MPI_Datatype MpiData<unsigned int long>::Type = MPI_LONG;

	template <> const MPI_Datatype MpiData<unsigned int>::Type = MPI_INTEGER;

	template <> const MPI_Datatype MpiData<double>::Type = MPI_DOUBLE;

	template <> const MPI_Datatype MpiData<int>::Type = MPI_INTEGER;

	void checkError(int errorCode, PsimagLite::String caller, CommType comm)
	{
		if (errorCode == MPI_SUCCESS)
			return;

		char errorMessage[MPI_MAX_ERROR_STRING];
		int  messageLength = 0;
		MPI_Error_string(errorCode, errorMessage, &messageLength);
		std::cerr << "Error in call to " << caller << " ";
		std::cerr << errorMessage << "\n";
		MPI_Abort(comm, -1);
	}

	bool hasMpi() { return true; }

	void init(int* argc, char** argv[]) { MPI_Init(argc, argv); }

	void finalize() { MPI_Finalize(); }

	SizeType commSize(CommType mpiComm)
	{
		int tmp = 1;
		if (mpiComm != 0)
			MPI_Comm_size(mpiComm, &tmp);
		return tmp;
	}

	SizeType commRank(CommType mpiComm)
	{
		int tmp = 0;
		if (mpiComm != 0)
			MPI_Comm_rank(mpiComm, &tmp);
		return tmp;
	}

	int barrier(CommType comm) { return MPI_Barrier(comm); }

	void version(std::ostream& os)
	{
		int version    = 0;
		int subversion = 0;
		int ret        = MPI_Get_version(&version, &subversion);
		checkError(ret, "MPI_Get_version");

		os << "MPI version=" << version << "." << subversion << "\n";

		int   resultlen  = 0;
		char* versionstr = new char[MPI_MAX_LIBRARY_VERSION_STRING + 1];
		ret              = MPI_Get_library_version(versionstr, &resultlen);
		checkError(ret, "MPI_Get_library_version");
		os << "MPI_Get_library_version=" << versionstr << "\n";

		delete[] versionstr;
		versionstr = nullptr;

		char* name = new char[MPI_MAX_PROCESSOR_NAME + 1];
		ret        = MPI_Get_processor_name(name, &resultlen);
		checkError(ret, "MPI_Get_processor_name");
		os << "MPI_Get_processor_name=" << name << "\n";

		delete[] name;
		name = nullptr;
	}

	void info(std::ostream& os)
	{
		int nkeys = 0;
		int ret   = MPI_Info_get_nkeys(MPI_INFO_ENV, &nkeys);
		checkError(ret, "MPI_Info_get_nkeys");

		if (nkeys == 0)
			return;

		char* key   = new char[MPI_MAX_INFO_KEY + 1];
		char* value = new char[MPI_MAX_INFO_VAL + 1];

		for (int i = 0; i < nkeys; ++i) {
			ret = MPI_Info_get_nthkey(MPI_INFO_ENV, i, key);
			checkError(ret, "MPI_Info_get_nthkey");

			int valuelen = 0;
			int flag     = 0;
			ret          = MPI_Info_get_valuelen(MPI_INFO_ENV, key, &valuelen, &flag);
			checkError(ret, "MPI_Info_get_valuelen");

			if (!flag || valuelen <= 0 || valuelen >= MPI_MAX_INFO_VAL)
				continue;

			ret = MPI_Info_get(MPI_INFO_ENV, key, valuelen, value, &flag);
			checkError(ret, "MPI_Info_Get");

			os << "MPI_INFO_ENV key=" << key << " value=" << value << "\n";
		}

		delete[] key;
		key = nullptr;
		delete[] value;
		value = nullptr;
	}

} // namespace MPI

} // namespace PsimagLite
