#pragma once
#include <cstddef>
#include <utility>

const std::size_t DAC_GRANULAR = 300;

class Matrix{
    public:
    Matrix(int,int);
    Matrix(int n, int m, double* d) : N(n), M(m), data(d){}
    ~Matrix();
    std::pair<int,int> shape();
    const Matrix mult(const Matrix&);
    const Matrix add(const Matrix&);
    const Matrix sub(const Matrix&);
    const Matrix T();
    double& at(int,int);
    private:
    std::size_t N;
    std::size_t M;
    double* data;
    const void multDAC(const double*, double*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t);
    const void multSeq(const double*, double*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t);

    const void TDAC(double*, std::size_t, std::size_t, std::size_t, std::size_t);
    const void TSeq(double*, std::size_t, std::size_t, std::size_t, std::size_t);
    const void _add(const double*, double*, std::size_t, int);
};