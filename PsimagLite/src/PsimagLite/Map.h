// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

/** \ingroup PsimagLite */
/*@{*/

/*!
 *
 *
 */
#ifndef MAP_HEADER_H
#define MAP_HEADER_H
#include "AllocatorCpu.h"
#include <map>

namespace PsimagLite {

template <typename Key, typename T, typename Compare = std::less<Key>> class Map {
public:

	using Type = std::map<Key, T, Compare, typename Allocator<std::pair<const Key, T>>::Type>;
}; // class Map

template <typename T> class IsMapLike {
public:

	enum
	{
		True = false
	};
};

template <typename Key, typename T, typename Compare>
class IsMapLike<std::map<Key, T, Compare, typename Allocator<std::pair<const Key, T>>::Type>> {
public:

	enum
	{
		True = true
	};
};

template <typename MapType>
typename EnableIf<IsMapLike<MapType>::True, void>::Type
printMap(std::ostream& os, const MapType& x, const String& label)
{
	using MapIteratorType = typename MapType::const_iterator;
	for (MapIteratorType it = x.begin(); it != x.end(); ++it) {
		os << label << "[" << it->first << "]=" << it->second << "\n";
	}
}

} // namespace PsimagLite

/*@}*/
#endif // MAP_HEADER_H
