// SPDX-FileCopyrightText: Copyright (c) 2009-2016, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]


/*! \file ProgressIndicator.cpp
 *
 *  This class handles output to a progress indicator (usually the terminal)
 */

#include "ProgressIndicator.h"

namespace PsimagLite {

MemoryUsage   ProgressIndicator::musage_;
OstringStream ProgressIndicator::buffer_(std::cout.precision());
bool          ProgressIndicator::bufferActive_ = false;

} // namespace PsimagLite
