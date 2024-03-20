#pragma once
#include <cstddef>
#include <utility>


const std::size_t PAR_GRANULAR = 10;
const std::size_t DAC_GRANULAR = 100;

class Matrix{
    public:
    Matrix(int,int);
    ~Matrix();
    std::pair<int,int> shape();
    void reshape(int,int);
    const Matrix mult(const Matrix&);
    // Matrix add(Matrix);
    // Matrix sub(Matrix);
    // Matrix Transpose();
    double& at(int,int);
    private:
    std::size_t N;
    std::size_t M;
    double* data;
    const void DACMult(const double*, double*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t);
    const void multSeq(const double*, double*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t);
    const void multPar(const double*, double*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t);
};