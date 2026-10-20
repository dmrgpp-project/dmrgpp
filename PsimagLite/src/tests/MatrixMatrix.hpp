#ifndef PSIMAGLITE_TEST_MATRIX_MATRIX_HPP
#define PSIMAGLITE_TEST_MATRIX_MATRIX_HPP

#include <PsimagLite/BLAS.h>

#include <Kokkos_Core.hpp>

#include <array>
#include <cassert>

class MatrixMatrix {
public:
    static constexpr SizeType Tasks = 4;

    static constexpr int M = 1024;
    static constexpr int N = 1024;
    static constexpr int K = 1024;

    using LeftMatrix   = std::vector<double>;
    using RightMatrix  = std::vector<double>;
    using OutputMatrix = std::vector<double>;

    MatrixMatrix()
    {
        for (SizeType task = 0; task < Tasks; ++task) {
            left_[task].resize(M * K);
            right_[task].resize(K * N);
            output_[task].resize(M * N);

            for (int column = 0; column < K; ++column) {
                for (int row = 0; row < M; ++row) {
                    left_[task][row + column * M] =
                        static_cast<double>(
                            100 * task + 10 * column + row + 1);
                }
            }

            for (int column = 0; column < N; ++column) {
                for (int row = 0; row < K; ++row) {
                    right_[task][row + column * K] =
                        static_cast<double>(
                            50 * task + 8 * column + row + 1);
                }
            }

            output_[task].assign(M * N, -999.0);
        }

        std::vector<int> weights(Tasks, 1);
        execs = Kokkos::Experimental::partition_space(
            Kokkos::DefaultExecutionSpace{}, weights);
    }

    SizeType tasks() const { return Tasks; }

    void doTask(SizeType taskNumber, SizeType)
    {
        std::cout << "Start task " << taskNumber << std::endl;
        assert(taskNumber < Tasks);

        struct ActiveGuard {
            std::atomic<int>& active;
            std::atomic<int>& max_active;

            ActiveGuard(std::atomic<int>& active_,
                        std::atomic<int>& max_active_)
                : active(active_)
                , max_active(max_active_)
            {
                const int current = active.fetch_add(1) + 1;

                int previous = max_active.load();
                while (previous < current &&
                       !max_active.compare_exchange_weak(
                           previous, current)) {
                }
            }

            ~ActiveGuard()
            {
                active.fetch_sub(1);
            }
        };

        ActiveGuard guard(active, max_active);

        PsimagLite::kokkos_gemm(
            execs[taskNumber],
            'N', 'N',
            M, N, K,
            1.0,
            left_[taskNumber].data(), M,
            right_[taskNumber].data(), K,
            0.0,
            output_[taskNumber].data(), M);
        std::cout << "End task " << taskNumber << std::endl;
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

    ~MatrixMatrix()
    {
        std::cout << "max_active: " << max_active.load()
                  << '\n';
    }

private:
    std::array<LeftMatrix, Tasks> left_;
    std::array<RightMatrix, Tasks> right_;
    std::array<OutputMatrix, Tasks> output_;

    std::atomic<int> active{0};
    std::atomic<int> max_active{0};

    // Keep this declaration compatible with the actual return type
    // of partition_space in the Kokkos version being tested.
    decltype(Kokkos::Experimental::partition_space(
        Kokkos::DefaultExecutionSpace{},
        std::vector<int>(Tasks, 1))) execs;
};

#endif // PSIMAGLITE_TEST_MATRIX_MATRIX_HPP
