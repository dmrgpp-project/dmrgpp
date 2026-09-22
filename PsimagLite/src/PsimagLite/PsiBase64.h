/*
 *****
 ***** ATTENTION: This code has been modified from the original *****
 *****            and adapted for use in PsimagLite
 *****
------------------------- begin original notice -------------------------------

Copyright (C) 2004-2008 René Nyffenegger

This source code is provided 'as-is', without any express or implied
warranty. In no event will the author be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this source code must not be misrepresented; you must not
   claim that you wrote the original source code. If you use this source code
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original source code.

3. This notice may not be removed or altered from any source distribution.

René Nyffenegger rene.nyffenegger@adp-gmbh.ch

------------------------- end original notice -------------------------------

***** Original notice is above
***** ATTENTION: This code has been modified from the original *****
*/

#ifndef PSIBASE64_H
#define PSIBASE64_H

#include "Io/IoSerializerStub.h"
#include "Vector.h"

namespace PsimagLite {

class PsiBase64 {

	static const String BASE64_CHARS;

public:

	class Encode {

	public:

		Encode(const String& str)
		{
			encode(reinterpret_cast<const unsigned char*>(str.c_str()), str.length());
		}

		Encode(unsigned char const* bytesToEncode, unsigned int inLen)
		{
			encode(bytesToEncode, inLen);
		}

		const String& operator()() const { return buffer_; }

		void write(String label, PsimagLite::IoSerializer& ioSerializer) const
		{
			ioSerializer.write(label, buffer_);
		}

		friend std::ostream& operator<<(std::ostream& os, const Encode& encode)
		{
			os << "#InputStartsHere\n";
			os << encode.buffer_;
			os << "\n#InputEndsHere\n";
			return os;
		}

	private:

		void encode(unsigned char const* bytesToEncode, unsigned int inLen)
		{
			buffer_         = "";
			int           i = 0;
			int           j = 0;
			unsigned char char_array3[3];
			unsigned char char_array4[4];

			while (inLen--) {
				char_array3[i++] = *(bytesToEncode++);
				if (i == 3) {
					char_array4[0] = (char_array3[0] & 0xfc) >> 2;
					char_array4[1] = ((char_array3[0] & 0x03) << 4)
					    + ((char_array3[1] & 0xf0) >> 4);
					char_array4[2] = ((char_array3[1] & 0x0f) << 2)
					    + ((char_array3[2] & 0xc0) >> 6);
					char_array4[3] = char_array3[2] & 0x3f;

					for (i = 0; (i < 4); i++)
						buffer_ += BASE64_CHARS[char_array4[i]];
					i = 0;
				}
			}

			if (i) {
				for (j = i; j < 3; j++)
					char_array3[j] = '\0';

				char_array4[0] = (char_array3[0] & 0xfc) >> 2;
				char_array4[1] = ((char_array3[0] & 0x03) << 4)
				    + ((char_array3[1] & 0xf0) >> 4);
				char_array4[2] = ((char_array3[1] & 0x0f) << 2)
				    + ((char_array3[2] & 0xc0) >> 6);
				char_array4[3] = char_array3[2] & 0x3f;

				for (j = 0; (j < i + 1); j++)
					buffer_ += BASE64_CHARS[char_array4[j]];

				while ((i++ < 3))
					buffer_ += '=';
			}
		}

		String buffer_;
	}; // class Encode

	class Decode {

	public:

		Decode(const String& encodedString)
		{
			buffer_              = "";
			int           in_len = encodedString.size();
			int           i      = 0;
			int           j      = 0;
			int           in     = 0;
			unsigned char char_array4[4], char_array3[3];

			while (in_len-- && (encodedString[in] != '=')
			       && isBase64(encodedString[in])) {
				char_array4[i++] = encodedString[in];
				in++;
				if (i == 4) {
					for (i = 0; i < 4; i++)
						char_array4[i] = BASE64_CHARS.find(char_array4[i]);

					char_array3[0] = (char_array4[0] << 2)
					    + ((char_array4[1] & 0x30) >> 4);
					char_array3[1] = ((char_array4[1] & 0xf) << 4)
					    + ((char_array4[2] & 0x3c) >> 2);
					char_array3[2]
					    = ((char_array4[2] & 0x3) << 6) + char_array4[3];

					for (i = 0; (i < 3); i++)
						buffer_ += char_array3[i];
					i = 0;
				}
			}

			if (i) {
				for (j = i; j < 4; j++)
					char_array4[j] = 0;

				for (j = 0; j < 4; j++)
					char_array4[j] = BASE64_CHARS.find(char_array4[j]);

				char_array3[0]
				    = (char_array4[0] << 2) + ((char_array4[1] & 0x30) >> 4);
				char_array3[1] = ((char_array4[1] & 0xf) << 4)
				    + ((char_array4[2] & 0x3c) >> 2);
				char_array3[2] = ((char_array4[2] & 0x3) << 6) + char_array4[3];

				for (j = 0; (j < i - 1); j++)
					buffer_ += char_array3[j];
			}
		}

		const String& operator()() const { return buffer_; }

	private:

		static bool isBase64(unsigned char c)
		{
			return (isalnum(c) || (c == '+') || (c == '/'));
		}

		String buffer_;
	}; // class Decode
}; // class PsiBase64
} // namespace PsimagLite

#endif // PSIBASE64_H
