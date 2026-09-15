// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Fermi.h
 *
 *  Fermi functions and their derivatives
 */

#ifndef FERMI_H_
#define FERMI_H_
#include <cmath>

namespace PsimagLite {

template <typename FieldType> FieldType fermi(const FieldType& x)
{
	if (x > 50)
		return 0;
	if (x < -50)
		return 1;
	if (x < 0)
		return 1.0 / (1.0 + std::exp(x));
	return std::exp(-x) / (1.0 + std::exp(-x));
}

// Derivative (prime) of Fermi's function
template <typename FieldType> FieldType fermiPrime(const FieldType& x)
{
	FieldType res;
	res = -fermi(x) * fermi(-x);
	return res;
}

template <typename FieldType> FieldType logfermi(const FieldType& x)
{
	FieldType res;
	if (x > 20)
		return -x;
	if (x < -20)
		return 0;
	res = -log(1.0 + exp(x));
	return res;
}

} // namespace PsimagLite

/*@}*/
#endif // FERMI_H_
