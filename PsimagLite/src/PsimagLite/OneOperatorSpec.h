#ifndef ONEOPERATORSPEC_H
#define ONEOPERATORSPEC_H
#include "PsimagLite.h"
#include "Vector.h"
#include <cstdlib>

namespace PsimagLite {
struct OneOperatorSpec {

	OneOperatorSpec(PsimagLite::String label_)
	    : dof(0)
	    , label(label_)
	    , transpose(false)
	{
		SizeType last_index = label.length();
		if (last_index > 0)
			last_index--;
		if (label[last_index] == '\'') {
			label     = label.substr(0, last_index);
			transpose = true;
		}

		label_ = label;

		SizeType i = 0;
		for (; i < label.length(); ++i) {
			if (label[i] == '?')
				break;
		}

		if (i == label.length())
			return;

		if (i + 1 == label.length())
			err("WRONG op. spec. " + label_ + ", nothing after ?\n");

		label                       = label_.substr(0, i);
		const String numeric_string = label_.substr(i + 1, label_.length());
		if (!isAnInteger(numeric_string)) {
			throw RuntimeError("FATAL: Syntax Error: The label " + label
			                   + " must be followed by an integer " + "and not "
			                   + numeric_string + "\n");
		}

		dof = atoi(numeric_string.c_str());
	}

	struct SiteSplit {

		SiteSplit(bool hasSiteString_, String root_, String siteString_)
		    : hasSiteString(hasSiteString_)
		    , root(root_)
		    , siteString(siteString_)
		{ }

		bool   hasSiteString;
		String root;
		String siteString;
	};

	static SiteSplit
	extractSiteIfAny(PsimagLite::String name, const char cBegin = '[', const char cEnd = ']')
	{
		int first_index = -1;
		int last_index  = -1;
		for (SizeType i = 0; i < name.length(); ++i) {
			if (name[i] == cBegin) {
				first_index = i;
				continue;
			}

			if (name[i] == cEnd) {
				last_index = i;
				continue;
			}
		}

		if (first_index < 0 && last_index < 0)
			return SiteSplit(false, name, "");

		bool b1 = (first_index < 0 && last_index >= 0);
		bool b2 = (first_index >= 0 && last_index < 0);
		if (b1 || b2)
			err(name + " has unmatched " + cBegin + " or " + cEnd + "\n");

		String str = name.substr(0, first_index);
		str += name.substr(last_index + 1, name.length() - last_index);
		String site_string = name.substr(first_index + 1, last_index - first_index - 1);
		return SiteSplit(true, str, site_string);
	}

	static bool isNonNegativeInteger(const String& s)
	{
		return !s.empty()
		    && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c); })
		    == s.end();
	}

	static SizeType strToNumberOrFail(String s)
	{
		if (!isNonNegativeInteger(s))
			err("string " + s + " is not a NonNegativeInteger\n");
		return atoi(s.c_str());
	}

	SizeType dof;
	String   label;
	bool     transpose;
}; // struct OneOperatorSpec

} // namespace PsimagLite
#endif // ONEOPERATORSPEC_H
