// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/*! \file ProgramGlobals.cpp
 *
 *
 *
 */

#include "ProgramGlobals.h"
#include "../Version.h"

namespace Dmrg {

SizeType                 ProgramGlobals::maxElectronsOneSpin = 0;
bool                     ProgramGlobals::oldChangeOfBasis    = false;
const PsimagLite::String ProgramGlobals::license
    = "Copyright (c) 2009-2016-2018, UT-Battelle, LLC\n"
      "\n"
      "[DMRG++, Version "
    + ttos(DMRGPP_VERSION)
    + "]\n"
      "\n"
      "---------------------------------------------------------\n"
      "Licensed under GNU GPL version 3 or later (GPL-3.0-or-later).\n"
      "This program comes with ABSOLUTELY NO WARRANTY.\n"
      "\n"
      "See LICENSE, NOTICE, and THIRD_PARTY_NOTICES distributed with\n"
      "this software.\n"
      "---------------------------------------------------------\n"
      "\n";

PsimagLite::String ProgramGlobals::SYSTEM_STACK_STRING  = "SystemStack";
PsimagLite::String ProgramGlobals::ENVIRON_STACK_STRING = "EnvironStack";

} // namespace Dmrg
