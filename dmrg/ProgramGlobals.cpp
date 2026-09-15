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
      "All rights reserved\n"
      "\n"
      "[DMRG++, Version "
    + ttos(DMRGPP_VERSION)
    + "]\n"
      "\n"
      "---------------------------------------------------------\n"
      "THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
      "CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
      "WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
      "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
      "PARTICULAR PURPOSE ARE DISCLAIMED.\n"
      "\n"
      "Please see full open source license included in file LICENSE.\n"
      "---------------------------------------------------------\n"
      "\n";

PsimagLite::String ProgramGlobals::SYSTEM_STACK_STRING  = "SystemStack";
PsimagLite::String ProgramGlobals::ENVIRON_STACK_STRING = "EnvironStack";

} // namespace Dmrg
