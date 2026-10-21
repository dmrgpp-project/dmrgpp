#ifndef PSIMAGLITE_TEST_MATRIX_MATRIX_HPP
#define PSIMAGLITE_TEST_MATRIX_MATRIX_HPP

#include <PsimagLite/BLAS.h>
#include <array>
#include <cassert>

class MatrixMatrix {
public:

	static constexpr SizeType Tasks = 4;
	using LeftMatrix                = std::array<double, 3 * 8>;
	using RightMatrix               = std::array<double, 8 * 3>;
	using OutputMatrix              = std::array<double, 3 * 3>;

	MatrixMatrix()
	{
		for (SizeType task = 0; task < Tasks; ++task) {
			for (SizeType column = 0; column < 8; ++column) {
				for (SizeType row = 0; row < 3; ++row)
					left_[task][row + column * 3] = static_cast<double>(
					    100 * task + 10 * column + row + 1);
			}

			for (SizeType column = 0; column < 3; ++column) {
				for (SizeType row = 0; row < 8; ++row)
					right_[task][row + column * 8]
					    = static_cast<double>(50 * task + 8 * column + row + 1);
			}

			output_[task].fill(-999.0);
		}
	}

	SizeType tasks() const { return Tasks; }

	void doTask(SizeType taskNumber, SizeType)
	{
		assert(taskNumber < Tasks);
		psimag::BLAS::GEMM('N',
		                   'N',
		                   3,
		                   3,
		                   8,
		                   1.0,
		                   left_[taskNumber].data(),
		                   3,
		                   right_[taskNumber].data(),
		                   8,
		                   0.0,
		                   output_[taskNumber].data(),
		                   3);
	}

	const LeftMatrix& left(SizeType taskNumber) const
	{
		assert(taskNumber < Tasks);
		return left_[taskNumber];
	}

	const RightMatrix& right(SizeType taskNumber) const
	{
		assert(taskNumber < Tasks);
		return right_[taskNumber];
	}

	const OutputMatrix& output(SizeType taskNumber) const
	{
		assert(taskNumber < Tasks);
		return output_[taskNumber];
	}

private:

	std::array<LeftMatrix, Tasks>   left_;
	std::array<RightMatrix, Tasks>  right_;
	std::array<OutputMatrix, Tasks> output_;
};

#endif // PSIMAGLITE_TEST_MATRIX_MATRIX_HPP
