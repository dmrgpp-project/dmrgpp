#ifndef PSI_MEM_RESOLV_H
#define PSI_MEM_RESOLV_H
#include "IsClass.h"
#include "Map.h"
#include "Sort.h"
#include "Vector.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>

namespace PsimagLite {

struct MemoryPointer {
	unsigned int      type;
	unsigned int      length;
	long unsigned int ptr;
}; // struct MemoryPointers

bool operator==(const MemoryPointer& a, const MemoryPointer& b);

class MemResolv {

	using VectorMemoryPointerType = std::vector<MemoryPointer>;
	using PairType                = std::pair<long unsigned int, long int>;
	using VectorPairType          = std::vector<PairType>;
	typedef double (*RefFunctionType)(double);

public:

	enum MemoryKindEnum
	{
		MEMORY_DATA,
		MEMORY_HEAPPTR,
		MEMORY_TEXTPTR
	};

	static const unsigned int      sizeof_heapref = sizeof(void*);
	static const unsigned int      sizeof_vptr    = sizeof(void*);
	static const unsigned int      sizeof_heapptr = sizeof(void*);
	static const long unsigned int label_length   = 128;

	template <typename T>
	MemResolv(T* ptr)
	    : intoOffset_(0)
	    , refTextPtr_(0)
	    , zeroes_(0)
	    , lenOfZeroes_(0)
	{ }

	MemResolv(String filename, String label)
	    : intoOffset_(0)
	    , refTextPtr_(0)
	    , zeroes_(0)
	    , lenOfZeroes_(0)
	{
		if (label.length() < label_length) {
			SizeType to_add = label_length - label.length();
			for (SizeType i = 0; i < to_add; ++i)
				label.push_back(0);
		}

		String mresolv_name("MemResolv::ctor():");

		std::ifstream fin(filename.c_str());
		if (!fin || fin.bad() || !fin.good())
			throw RuntimeError(mresolv_name + " cannot open " + filename + "\n");

		long unsigned int len_of_label = 0;
		fin.read(reinterpret_cast<char*>(&len_of_label), sizeof(len_of_label));
		fin.read(reinterpret_cast<char*>(&len_of_label), sizeof(len_of_label));
		if (len_of_label != label_length)
			throw RuntimeError(mresolv_name + " label length error\n");
		if (len_of_label != label.length())
			throw RuntimeError(mresolv_name + " mismatched label length\n");

		String label2;
		label2.resize(label_length);
		fin.read(const_cast<char*>(label2.data()), label_length);
		if (!stringEqual(label2, label))
			throw RuntimeError(mresolv_name + " mismatched label");

		long unsigned int old_start     = 0;
		char*             ptr_old_start = reinterpret_cast<char*>(&old_start);
		fin.read(ptr_old_start, sizeof(old_start));

		std::cout << "Recovered reference heap pointer ";
		std::cout << reinterpret_cast<void*>(old_start) << "\n";

		fin.read(reinterpret_cast<char*>(&refTextPtr_), sizeof(refTextPtr_));
		std::cout << "Recovered reference text pointer " << refTextPtr_ << "\n";

		fin.read(reinterpret_cast<char*>(&intoOffset_), sizeof(intoOffset_));

		loadChunkInfo(fin);

		long unsigned int len     = 0;
		char*             ptr_len = reinterpret_cast<char*>(&len);
		fin.read(ptr_len, sizeof(len));
		std::cout << "MemResolv read from file len= " << len << "\n";
		unsigned char* source_ptr = new unsigned char[len];
		garbage_.push_back(source_ptr);
		garbageSize_.push_back(len);
		fin.read(reinterpret_cast<char*>(source_ptr), len);
		fin.close();

		// ADJUST POINTER VALUES
		long int new_start = pointerToLui(reinterpret_cast<void*>(source_ptr));
		long int offset    = new_start - old_start;

		adjustPointers(source_ptr, offset);
	}

	~MemResolv()
	{
		delete[] zeroes_;
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			unsigned char* ptr = garbage_[i];
			if (garbageSize_[i] > 0)
				delete[] ptr;
		}
	}

	void save(String filename, String label) const
	{
		assert(garbage_.size() > 0);

		if (label.length() < label_length) {
			SizeType to_add = label_length - label.length();
			for (SizeType i = 0; i < to_add; ++i)
				label.push_back(0);
		}

		if (label.length() != label_length)
			throw RuntimeError("MemResolv::save(): label length\n");

		std::ofstream fout(filename.c_str());

		if (!fout) {
			String msg("MemResolv::save(): cannot open file ");
			throw RuntimeError(msg + filename + "\n");
		}

		long unsigned int len_of_label = label.length();
		assert(lenOfLabel == LABEL_LENGTH);
		fout.write(reinterpret_cast<char*>(&len_of_label), sizeof(len_of_label));
		fout.write(reinterpret_cast<char*>(&len_of_label), sizeof(len_of_label));
		fout.write(label.data(), len_of_label);

		SizeType       total         = 0;
		SizeType       max_hole_size = 0;
		VectorPairType offsets_for_holes;
		findSizes(total, max_hole_size, offsets_for_holes);

		long unsigned int ref_ptr_value
		    = pointerToLui(reinterpret_cast<void*>(vmptr_[0].ptr));
		adjustPointer(reinterpret_cast<unsigned char*>(&ref_ptr_value),
		              sizeof(ref_ptr_value),
		              0,
		              &offsets_for_holes);
		char* ptr_ref_ptr = reinterpret_cast<char*>(&ref_ptr_value);
		fout.write(ptr_ref_ptr, sizeof(ref_ptr_value));

		fout.write(reinterpret_cast<const char*>(&refTextPtr_), sizeof(refTextPtr_));
		std::cout << "Written reference text pointer " << refTextPtr_ << "\n";

		const char* ptr_into_ptr = reinterpret_cast<const char*>(&intoOffset_);
		fout.write(ptr_into_ptr, sizeof(intoOffset_));

		updateZeroes(max_hole_size + 1);

		saveChunkInfo(fout, offsets_for_holes);

		long unsigned int len     = total;
		char*             ptr_len = reinterpret_cast<char*>(&len);
		fout.write(ptr_len, sizeof(len));

		SizeType total2 = saveChunkData(fout, offsets_for_holes);
		std::cout << "Saved " << total2 << " bytes to " << filename << "\n";
		fout.close();
	}

	unsigned char* get() { return garbage_[0] + intoOffset_; }

	template <typename T>
	void push(MemoryKindEnum type, unsigned int length, T* ptr, String msg = "")
	{
		MemoryPointer mptr;
		mptr.type        = type;
		mptr.length      = length;
		const void* vptr = reinterpret_cast<const void*>(ptr);
		mptr.ptr         = pointerToLui(vptr);
		if (find(vmptr_.begin(), vmptr_.end(), mptr) != vmptr_.end())
			return;

		vmptr_.push_back(mptr);
		rankVector_.push_back(mptr.ptr);

		if (msg == "")
			return;
		std::cout << "Address start " << vptr << " length " << mptr.length;
		std::cout << " type=" << type << " " << msg << "\n";
	}

	unsigned char* dup()
	{
		SizeType       total         = 0;
		SizeType       max_hole_size = 0;
		VectorPairType offsets_for_holes;
		findSizes(total, max_hole_size, offsets_for_holes);
		std::cout << "total = " << total << " maxHoleSize= " << max_hole_size << "\n";
		unsigned char* ptr = new unsigned char[total];
		garbage_.push_back(ptr);
		garbageSize_.push_back(total);

		updateZeroes(max_hole_size + 1);

		deepCopy(ptr, total, offsets_for_holes);

		return ptr + intoOffset_;
	}

	template <typename NoClassType>
	typename EnableIf<!IsClass<NoClassType>::value, SizeType>::Type
	memResolv(const NoClassType* c, SizeType x = sizeof(NoClassType), String msg = "")
	{
		SizeType tmp = sizeof(NoClassType);
		assert(x >= tmp);
		SizeType r     = x - tmp;
		char*    c_ptr = (char*)c;
		updateZeroes(r + 1, 0);
		c_ptr += tmp;
		memcpy(c_ptr, zeroes_, r);
		tmp += r;
		push(MemResolv::MEMORY_DATA, tmp, c, msg + " fundamental type");
		return tmp;
	}

	template <typename SomeComplexNumber>
	typename EnableIf<IsComplexNumber<SomeComplexNumber>::True, SizeType>::Type
	memResolv(const SomeComplexNumber*, SizeType = sizeof(SomeComplexNumber), String = "")
	{
		throw RuntimeError("memResolv for std::complex not implemented\n");
	}

	template <typename T1, typename T2, typename T3>
	SizeType memResolv(const std::basic_string<T1, T2, T3>* c,
	                   SizeType x   = sizeof(std::basic_string<T1, T2, T3>),
	                   String   msg = "")
	{
		assert(x == sizeof(std::basic_string<T1, T2, T3>));
		SizeType total = x;
		push(MemResolv::MEMORY_HEAPPTR, x, c, msg + " string ptr");

		const long unsigned int* c2  = (const long unsigned int*)c;
		const char*              c3  = (const char*)*c2;
		SizeType                 len = sizeof(unsigned int);
		char*                    ptr = const_cast<char*>(c3);
		ptr -= len;
		updateZeroes(len + 1, 0);
		memcpy(reinterpret_cast<char*>(ptr), zeroes_, len);
		ptr += len;
		len *= 3;
		ptr -= len;
		push(MemResolv::MEMORY_DATA, len, ptr, msg + " string misc");
		if (c->length() == 0)
			return total;
		push(MemResolv::MEMORY_DATA, c->length(), c3, msg + " string data");
		return total;
	}

	template <typename SomePairType>
	typename EnableIf<IsPairLike<SomePairType>::True, SizeType>::Type
	memResolv(const SomePairType* c, SizeType x = sizeof(SomePairType), String msg = "")
	{
		using FirstType  = typename SomePairType::first_type;
		using SecondType = typename SomePairType::second_type;

		const FirstType* f     = (const FirstType*)c;
		SizeType         total = memResolv(f, sizeof(FirstType), msg + "->first");

		const char* s = (const char*)c;
		s += sizeof(FirstType*);
		const SecondType* s2 = (const SecondType*)s;
		total += memResolv(s2, sizeof(SecondType), msg + "->second");

		assert(x >= sizeof(SomePairType));
		SizeType r = x - sizeof(SomePairType);

		if (r == 0)
			return total;

		updateZeroes(r + 1, 0);
		const char* s3  = reinterpret_cast<const char*>(c);
		char*       ptr = const_cast<char*>(s3);
		ptr += sizeof(SomePairType);
		memcpy(reinterpret_cast<char*>(ptr), zeroes_, r);

		return total;
	}

	template <typename SomeVectorType>
	typename EnableIf<IsVectorLike<SomeVectorType>::True
	                      && !IsClass<typename SomeVectorType::value_type>::value,
	                  SizeType>::Type
	memResolv(const SomeVectorType* v, SizeType = 0, String msg = "")
	{
		using SomeElementType = typename SomeVectorType::value_type;

		SizeType tmp = sizeof(SomeVectorType);
		push(MemResolv::MEMORY_HEAPPTR, tmp, v, msg + "vector<no class> ");
		assert(tmp = 3 * sizeof(void*)); // begin end reserve
		checkForVectorReserve(v);

		SizeType              total = tmp;
		const SomeVectorType& vv    = *v;

		tmp = sizeof(SomeElementType) * vv.size();
		if (tmp == 0)
			return total;

		push(MemResolv::MEMORY_DATA, tmp, &(vv[0]), msg + " vector data no class");
		// total += tmp;

		return total;
	}

	template <typename SomeVectorType>
	typename EnableIf<IsVectorLike<SomeVectorType>::True
	                      && IsClass<typename SomeVectorType::value_type>::value,
	                  SizeType>::Type
	memResolv(const SomeVectorType* v, SizeType = 0, String msg = "")
	{
		SizeType tmp = sizeof(SomeVectorType);
		assert(tmp = 3 * sizeof(void*)); // begin end reserve
		push(MemResolv::MEMORY_HEAPPTR, tmp, v, msg + " vector<class>");
		checkForVectorReserve(v);

		SizeType              total = tmp;
		const SomeVectorType& vv    = *v;

		if (vv.size() == 0)
			return total;
		SizeType element_size = sizeof(vv[0]);
		for (SizeType i = 0; i < vv.size(); ++i)
			memResolv(&vv[i], element_size, msg);

		return total;
	}

	template <typename SomeVectorType>
	typename EnableIf<IsVectorLike<SomeVectorType>::True
	                      && !IsClass<typename SomeVectorType::value_type>::value,
	                  SizeType>::Type
	memResolvPtr(const SomeVectorType* v, SizeType = 0, String msg = "")
	{
		using SomeElementType = typename SomeVectorType::value_type;

		SizeType tmp = sizeof(SomeVectorType);
		push(MemResolv::MEMORY_HEAPPTR, tmp, v, msg + " vector<no class> ");
		assert(tmp = 3 * sizeof(void*)); // begin end reserve
		checkForVectorReserve(v);

		SizeType              total = tmp;
		const SomeVectorType& vv    = *v;

		tmp = sizeof(SomeElementType) * vv.size();
		if (tmp == 0)
			return total;

		push(MemResolv::MEMORY_HEAPPTR, tmp, &(vv[0]), msg + " vector data is pointers");
		total += tmp;

		return total;
	}

	template <typename SomeMapType>
	typename EnableIf<IsMapLike<SomeMapType>::True, SizeType>::Type
	memResolv(const SomeMapType*, SizeType = 0, String = "")
	{
		throw RuntimeError("memResolv for std::map not implemented\n");
	}

	template <typename SomeClassType>
	typename EnableIf<!IsVectorLike<SomeClassType>::True && !IsPairLike<SomeClassType>::True
	                      && !IsMapLike<SomeClassType>::True
	                      && !IsComplexNumber<SomeClassType>::True
	                      && IsClass<SomeClassType>::value,
	                  SizeType>::Type
	memResolv(const SomeClassType* c, SizeType x = 0, String msg = "")
	{
		return c->memResolv(*this, x, msg);
	}

	friend std::ostream& operator<<(std::ostream& os, const MemResolv& mresolv);

private:

	template <typename SomeVectorType>
	void checkForVectorReserve(const SomeVectorType* vPtr) const
	{
		const long unsigned int* v       = (const long unsigned int*)vPtr;
		const long unsigned int* end     = (const long unsigned int*)(v + 1);
		const long unsigned int* reserve = (const long unsigned int*)(v + 2);
		if (*end != *reserve)
			std::cerr << "WARNING: std::vector has reserve\n";
	}

	void updateZeroes(SizeType x, int value = -128) const
	{
		if (lenOfZeroes_ == x)
			return;

		if (zeroes_ != 0)
			delete[] zeroes_;

		assert(x > 0);
		zeroes_ = new char[x];
		for (SizeType i = 0; i < x; ++i)
			zeroes_[i] = value;
		lenOfZeroes_ = x;
	}

	bool stringEqual(String str1, String str2) const
	{
		if (str1.length() != str2.length())
			return false;

		for (SizeType i = 0; i < str1.length(); ++i) {
			if (str1[i] == 0 && str2[i] == 0)
				break;
			if (str1[i] != str2[i])
				return false;
		}

		return true;
	}

	void print(std::ostream& os, const MemoryPointer& mp) const
	{
		os << mp.type << " " << mp.length << " " << mp.ptr << " 0x";
		std::hex(os);
		os << mp.ptr << "\n";
		std::dec(os);
	}

	void finish()
	{
		if (vmptr_.size() == 0)
			return;

		std::vector<SizeType>       iperm(rankVector_.size());
		Sort<std::vector<SizeType>> sort;
		sort.sort(rankVector_, iperm);
		unsigned int long old_start = pointerToLui(reinterpret_cast<void*>(vmptr_[0].ptr));
		VectorMemoryPointerType vmptr(vmptr_.size());
		std::cout << "Chunks in order\n";
		for (SizeType i = 0; i < iperm.size(); ++i) {
			SizeType j = iperm[i];
			vmptr[i]   = vmptr_[j];
			std::cout << "Entry " << vmptr[i].type << " 0x";
			std::hex(std::cout);
			std::cout << vmptr[i].ptr << " ";
			std::dec(std::cout);
			std::cout << " " << vmptr[i].length << "\n";
		}

		vmptr_                      = vmptr;
		unsigned int long new_start = pointerToLui(reinterpret_cast<void*>(vmptr_[0].ptr));
		intoOffset_
		    = (new_start > old_start) ? new_start - old_start : old_start - new_start;

		SizeType       total         = 0;
		SizeType       max_hole_size = 0;
		VectorPairType offsets_for_holes;
		findSizes(total, max_hole_size, offsets_for_holes);

		int long corrected_into_offset = intoOffset_;
		int long corrected_old_start   = old_start;
		adjustPointer(reinterpret_cast<unsigned char*>(&corrected_old_start),
		              sizeof(corrected_old_start),
		              0,
		              &offsets_for_holes);
		corrected_into_offset -= old_start;
		corrected_into_offset += corrected_old_start;

		intoOffset_ = corrected_into_offset;
	}

	void saveChunkInfo(std::ofstream& fout, const VectorPairType& offsetsForHoles) const
	{
		long unsigned int len = vmptr_.size();
		fout.write(reinterpret_cast<char*>(&len), sizeof(len));
		for (SizeType i = 0; i < vmptr_.size(); ++i) {
			MemoryPointer mptr = vmptr_[i];
			adjustPointer(reinterpret_cast<unsigned char*>(&(mptr.ptr)),
			              sizeof(mptr.ptr),
			              0,
			              &offsetsForHoles);

			fout.write(reinterpret_cast<char*>(&mptr), sizeof(MemoryPointer));
		}
	}

	SizeType saveChunkData(std::ofstream& fout, const VectorPairType& offsetsForHoles) const
	{
		SizeType total = 0;

		for (SizeType i = 0; i < vmptr_.size(); ++i) {
			if (i > 0) {
				long unsigned int end = vmptr_[i - 1].ptr + vmptr_[i - 1].length;
				SizeType          src_hole_size  = vmptr_[i].ptr - end;
				SizeType          dest_hole_size = src_hole_size % 8;

				if (lenOfZeroes_ <= dest_hole_size)
					throw RuntimeError("lenZeroes\n");
				if (dest_hole_size > 0) {
					fout.write(zeroes_, dest_hole_size);
				}

				total += dest_hole_size;
			}

			SizeType len       = vmptr_[i].length;
			char*    mptr      = reinterpret_cast<char*>(vmptr_[i].ptr);
			char*    allocated = 0;

			if (vmptr_[i].type == MEMORY_HEAPPTR) {
				allocated = new char[len];
				memcpy(allocated, reinterpret_cast<char*>(vmptr_[i].ptr), len);
				adjustPointer(reinterpret_cast<unsigned char*>(allocated),
				              len,
				              0,
				              &offsetsForHoles);
			}

			fout.write((allocated == 0) ? mptr : allocated, len);
			total += len;
			if (allocated)
				delete[] allocated;
		}

		return total;
	}

	void loadChunkInfo(std::ifstream& fin)
	{
		long unsigned int len = 0;
		fin.read(reinterpret_cast<char*>(&len), sizeof(len));
		vmptr_.resize(len);
		for (SizeType i = 0; i < vmptr_.size(); ++i) {
			MemoryPointer* mptr = &(vmptr_[i]);
			fin.read(reinterpret_cast<char*>(mptr), sizeof(MemoryPointer));
			rankVector_.push_back(mptr->ptr);
		}
	}

	void
	findSizes(SizeType& total, SizeType& maxHoleSize, VectorPairType& offsetsForHoles) const
	{
		maxHoleSize = 1;
		total       = 0;
		if (vmptr_.size() == 0)
			return;
		long unsigned int start = vmptr_[0].ptr;
		long unsigned int end   = start + vmptr_[0].length;
		total                   = vmptr_[0].length;
		for (SizeType i = 1; i < vmptr_.size(); ++i) {
			long unsigned int start2 = vmptr_[i].ptr;
			long unsigned int end2   = start2 + vmptr_[i].length;
			if (start2 < end)
				throw RuntimeError("findTotal end\n");

			total += vmptr_[i].length;
			long int src_hole_size  = start2 - end;
			long int dest_hole_size = src_hole_size % 8;
			long int offset         = dest_hole_size - src_hole_size;
			PairType offset_for_hole(start2, offset);
			offsetsForHoles.push_back(offset_for_hole);
			total += dest_hole_size;
			if (static_cast<long unsigned int>(dest_hole_size) > maxHoleSize)
				maxHoleSize = dest_hole_size;

			end   = end2;
			start = start2;
		}
	}

	void deepCopy(unsigned char* ptr, SizeType total, const VectorPairType& offsetsForHoles)
	{
		if (vmptr_.size() == 0)
			return;

		long unsigned int new_start = pointerToLui(reinterpret_cast<void*>(ptr));
		long unsigned int old_start = vmptr_[0].ptr;
		long int          offset    = new_start - old_start;

		long unsigned int start = vmptr_[0].ptr;
		SizeType          len   = vmptr_[0].length;
		long unsigned int end   = start + len;

		SizeType total2 = copyData(&ptr, 0, offset, offsetsForHoles);

		for (SizeType i = 1; i < vmptr_.size(); ++i) {
			if (total2 >= total)
				throw RuntimeError("deepCopy total2\n");
			long unsigned int start2 = vmptr_[i].ptr;
			SizeType          len    = vmptr_[i].length;
			long unsigned int end2   = start2 + len;
			if (start2 < end)
				throw RuntimeError("deepCopy\n");
			SizeType src_hole_size  = start2 - end;
			SizeType dest_hole_size = src_hole_size % 8;

			if (lenOfZeroes_ <= dest_hole_size)
				throw RuntimeError("lenZeroes\n");
			if (dest_hole_size > 0) {
				memcpy(ptr, zeroes_, dest_hole_size);
			}

			ptr += dest_hole_size;
			total2 += dest_hole_size;

			total2 += copyData(&ptr, i, offset, offsetsForHoles);

			start = start2;
			end   = end2;
		}

		if (total2 < total) {
			assert(total - total2 < lenOfZeroes_);
			memcpy(ptr, zeroes_, total - total2);
		}
	}

	SizeType copyData(unsigned char**       ptr,
	                  SizeType              i,
	                  long int              offset,
	                  const VectorPairType& offsetsForHoles)
	{
		const void* src  = reinterpret_cast<const void*>(vmptr_[i].ptr);
		void*       src2 = const_cast<void*>(src);
		SizeType    len  = vmptr_[i].length;

		memcpy(*ptr, src2, len);
		if (vmptr_[i].type == MEMORY_HEAPPTR)
			adjustPointer(*ptr, len, offset, &offsetsForHoles);
		*ptr += len;
		return len;
	}

	void adjustPointers(unsigned char*, long int) const
	{
		/*VectorPairType offsetsForHoles(0);
		RefFunctionType f = &PsimagLite::conj;
		long int offsetText = (long int)(*f);
		std::cout<<offsetText<<"\n";
		offsetText -= refTextPtr_;

		for (SizeType i = 0; i < vmptr_.size(); ++i) {
		        if (i > 0) {
		                long unsigned int end = vmptr_[i-1].ptr +
		vmptr_[i-1].length; SizeType holeSize = vmptr_[i].ptr - end; if
		(holeSize > 0) std::cerr<<"WARNING: input file has holes\n"; ptr
		+= holeSize;
		        }

		        SizeType len = vmptr_[i].length;
		        if (vmptr_[i].type == MEMORY_HEAPPTR)
		                adjustPointer(ptr,len,offset,&offsetsForHoles);
		        if (vmptr_[i].type == MEMORY_TEXTPTR)
		                adjustPointer(ptr,len,offsetText,0);
		        ptr += len;

		}*/
	}

	void adjustPointer(unsigned char*        ptr,
	                   SizeType              n,
	                   long int              offset,
	                   const VectorPairType* offsetsForHoles) const
	{
		if (n < 8 || n % 8 != 0)
			throw RuntimeError("adjustPointer");

		SizeType counter = 0;

		do {
			void*              p          = reinterpret_cast<void*>(ptr);
			long unsigned int* ptr_to_lui = reinterpret_cast<long unsigned int*>(ptr);
			long unsigned int  value      = *ptr_to_lui;

			if (value != 0) {
				long int correct_for_holes
				    = correctionForHoles(value, offsetsForHoles);
				long int all_offset_corrections = correct_for_holes + offset;
				value += all_offset_corrections;
				long unsigned int* value_ptr = &value;
				memcpy(p, value_ptr, 8);
			}

			ptr += 8;
			counter += 8;
		} while (counter < n);
	}

	long int correctionForHoles(long unsigned int     value,
	                            const VectorPairType* offsetsForHolesPtr) const
	{
		if (offsetsForHolesPtr == 0)
			return 0;
		const VectorPairType& offsets_for_holes = *offsetsForHolesPtr;
		long int              c                 = 0;
		for (SizeType i = 0; i < offsets_for_holes.size(); ++i) {
			SizeType start = offsets_for_holes[i].first;
			if (value < start)
				return c;
			c += offsets_for_holes[i].second;
		}

		return c;
	}

	long unsigned int pointerToLui(const void* ptr) const
	{
		return reinterpret_cast<long unsigned int>(ptr);
	}

	long unsigned int           intoOffset_;
	long int                    refTextPtr_;
	mutable char*               zeroes_;
	mutable SizeType            lenOfZeroes_;
	VectorMemoryPointerType     vmptr_;
	std::vector<SizeType>       rankVector_;
	std::vector<unsigned char*> garbage_;
	std::vector<SizeType>       garbageSize_;
}; // class MemResolv

std::ostream& operator<<(std::ostream& os, const MemResolv& mresolv);

template <typename T, bool IsClassTrait> class ResolveFinalOrNot { };

template <typename T> class ResolveFinalOrNot<T, false> {

public:

	static SizeType resolve(MemResolv& vmptr, SizeType size, const T* ptr)
	{
		SizeType tmp = sizeof(T) * size;
		vmptr.push(MemResolv::MEMORY_DATA, tmp, ptr, "pod");
		return tmp;
	}
};

template <typename T> class ResolveFinalOrNot<T, true> {

public:

	static SizeType resolve(MemResolv& vmptr, SizeType size, const T* ptr)
	{
		SizeType tmp = 0;
		for (SizeType i = 0; i < size; ++i)
			tmp += vmptr.memResolv(&(ptr[i]));
		return tmp;
	}
};

} // namespace PsimagLite

#endif // PSI_MEM_RESOLV_H
