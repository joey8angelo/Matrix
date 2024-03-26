#pragma once
#include <cstddef>
#include <utility>
#include <string>
#include "parallel.h"

const std::size_t DAC_GRANULAR = 300;

template<typename P>
class Matrix{
    public:
    Matrix(int,int);
    Matrix(int,int,P);
    Matrix(int n, int m, P* d) : N(n), M(m), data(d){} // construct with given data
    ~Matrix();
    std::pair<int,int> shape() const;
    Matrix<P> dot(const Matrix<P>&) const;
    Matrix<P> add(const Matrix<P>&) const;
    Matrix<P> sub(const Matrix<P>&) const;
    Matrix<P> mult(const Matrix<P>&) const;
    Matrix<P> T() const;
    P& at(int,int) const;
    template<typename U> Matrix<P> add(const U&) const;
    template<typename U> Matrix<P> sub(const U&) const;
    template<typename U> Matrix<P> mult(const U&) const;
    std::string strShape() const;
    void operator=(const Matrix<P>&);
    private:
    std::size_t N;
    std::size_t M;
    P* data;
    void dotDAC(const P*, P*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t) const;
    void dotSeq(const P*, P*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t) const;
    void TDAC(P*, std::size_t, std::size_t, std::size_t, std::size_t) const;
    void TSeq(P*, std::size_t, std::size_t, std::size_t, std::size_t) const;
};

/* Assignment operator */
template<typename P>
void Matrix<P>::operator=(const Matrix<P>& o){
    // if number of elements differ, reallocate memory
    if(N*M != o.N*o.M){
        delete[] data;
        data = new P[N * M];
    }
    // update size
    N = o.N;
    M = o.M;
    // populate data
    for(std::size_t i = 0; i < N * M; i++){
        data[i] = o.data[i];
    }
}
/* Element-wise addition of matrices */
template<typename P>
inline Matrix<P> operator+(const Matrix<P>& lhs, const Matrix<P>& rhs){
    return lhs.add(rhs);
}
/* Element-wise addition of a scalar on the left hand side */
template<typename P, typename U>
inline Matrix<P> operator+(const U& lhs, const Matrix<P>& rhs){
    return rhs.add(lhs);
}
/* Element-wise addition of a scalar on the right hand side */
template<typename P, typename U>
inline Matrix<P> operator+(const Matrix<P>& lhs, const U& rhs){
    return lhs.add(rhs);
}
/* Element-wise subtraction of matrices */
template<typename P>
inline Matrix<P> operator-(const Matrix<P>& lhs, const Matrix<P>& rhs){
    return lhs.sub(rhs);
}
/* Element-wise subtraction of a scalar on the left hand side */
template<typename P, typename U>
inline Matrix<P> operator-(const U& lhs, const Matrix<P>& rhs){
    return rhs.mult(-1).add(lhs);
}
/* Element-wise subtraction of a scalar on the right hand side */
template<typename P, typename U>
inline Matrix<P> operator-(const Matrix<P>& lhs, const U& rhs){
    return lhs.sub(rhs);
}
/* Element-wise multiplication of matrices */
template<typename P>
inline Matrix<P> operator*(const Matrix<P>& lhs, const Matrix<P>& rhs){
    return lhs.mult(rhs);
}
/* Element-wise multiplication of a scalar on the left hand side */
template<typename P, typename U>
inline Matrix<P> operator*(const U& lhs, const Matrix<P>& rhs){
    return rhs.mult(lhs);
}
/* Element-wise multiplication of a scalar on the right hand side */
template<typename P, typename U>
inline Matrix<P> operator*(const Matrix<P>& lhs, const U& rhs){
    return lhs.mult(rhs);
}

/* Default value Constructor */
template<typename P>
Matrix<P>::Matrix(int n, int m) : N(n), M(m){
    if(N <= 0 || M <= 0) throw std::runtime_error("Invalid matrix dimensions: " + strShape());
    data = new P[N * M];
    for(std::size_t i = 0; i < N * M; i++){
        data[i] = P();
    }
}

/* Constructor */
template<typename P>
Matrix<P>::Matrix(int n, int m, P v) : N(n), M(m){
    if(N <= 0 || M <= 0) throw std::runtime_error("Invalid matrix dimensions: " + strShape());
    data = new P[N * M];
    for(std::size_t i = 0; i < N * M; i++){
        data[i] = v;
    }
}

/* Destructor */
template<typename P>
Matrix<P>::~Matrix(){
    delete[] data;
}

/* Returns a pair with the number of rows and the number of columns in the matrix */
template<typename P>
std::pair<int,int> Matrix<P>::shape() const {
    return std::make_pair(N,M);
}

/* Returns a string representation of the shape of the matrix: "(N,M)" */
template<typename P>
std::string Matrix<P>::strShape() const {
    return "(" + std::to_string(N) + "," + std::to_string(M) + ")";
}

/* Returns a reference to the position in the matrix at the given coordinates */
template<typename P>
P& Matrix<P>::at(int i, int j) const {
    if(i < 0 || i >= N || j < 0 || j >= M) throw std::runtime_error("Index out of bounds: " + std::to_string(i) + " " + std::to_string(j));
    return data[i * M + j];
}

/* Computes the dot product of two matrices, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::dot(const Matrix<P>& B) const {
    if(M != B.N) throw std::runtime_error("Incompatible matrix dimensions for dot product: " + strShape() + " " + B.strShape());
    Matrix C(N, B.M);
    dotDAC(B.data, C.data, B.M, N, M, B.N, B.M, 0, 0, 0, 0, 0, 0);
    return C;
}

/* Element-wise addition of two matrices, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::add(const Matrix<P>& B) const {
    if(N != B.N || M != B.M) throw std::runtime_error("Incompatible matrix dimensions for addition: " + strShape() + " " + B.strShape());
    P* Cdata = new P[N * M];
    for(std::size_t i = 0; i < N*M; i++) Cdata[i] = data[i] + B.data[i];
    return Matrix<P>(N,M,Cdata);
}
/* Element-wise addition of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::add(const U& v) const {
    P* Cdata = new P[N * M];
    for(std::size_t i = 0; i < N*M; i++) Cdata[i] = data[i] + v;
    return Matrix<P>(N,M,Cdata);
}
/* Element-wise multiplication of two matrices, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::mult(const Matrix<P>& B) const {
    if(N != B.N || M != B.M) throw std::runtime_error("Incompatible matrix dimensions for multiplication: " + strShape() + " " + B.strShape());
    P* Cdata = new P[N * B.M];
    for(std::size_t i = 0; i < N*M; i++) Cdata[i] = data[i] * B.data[i];
    return Matrix<P>(N,B.M,Cdata);
}
/* Element-wise multiplication of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::mult(const U& v) const {
    P* Cdata = new P[N * M];
    for(std::size_t i = 0; i < N*M; i++) Cdata[i] = data[i] * v;
    return Matrix<P>(N,M,Cdata);
}
/* Element-wise subtraction of two matrices, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::sub(const Matrix<P>& B) const {
    if(N != B.N || M != B.M) throw std::runtime_error("Incompatible matrix dimensions for subtraction: " + strShape() + " " + B.strShape());
    P* Cdata = new P[N * M];
    for(std::size_t i = 0; i < N*M; i++) Cdata[i] = data[i] - B.data[i];
    return Matrix<P>(N,M,Cdata);
}
/* Element-wise subtraction of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::sub(const U& v) const {
    P* Cdata = new P[N * M];
    for(std::size_t i = 0; i < N*M; i++) Cdata[i] = data[i] - v;
    return Matrix<P>(N,M,Cdata);
}

/* Divide and Conquer inner product, allows for parallel execution and better IO efficiency than the basic inner product operation*/
template<typename P>
void Matrix<P>::dotDAC(const P* B, P* C, std::size_t BM, 
             std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, 
             std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM) const {
    if(std::min(Na, std::min(Ma, std::min(Nb, Mb))) <= DAC_GRANULAR) {
        dotSeq(B,C,BM,Na,Ma,Nb,Mb,PaN,PaM,PbN,PbM,PcN,PcM);
    }
    else {
        // size of submatrices
        std::size_t Nad = Na / 2;
        std::size_t Mad = Ma / 2;
        std::size_t Nbd = Nb / 2;
        std::size_t Mbd = Mb / 2;

        // positions of submatrices
        std::size_t A11N = PaN;
        std::size_t A11M = PaM;
        std::size_t A12N = PaN;
        std::size_t A12M = PaM + Mad;
        std::size_t A21N = PaN + Nad;
        std::size_t A21M = PaM;
        std::size_t A22N = PaN + Nad;
        std::size_t A22M = PaM + Mad;

        std::size_t B11N = PbN;
        std::size_t B11M = PbM;
        std::size_t B12N = PbN;
        std::size_t B12M = PbM + Mbd;
        std::size_t B21N = PbN + Nbd;
        std::size_t B21M = PbM;
        std::size_t B22N = PbN + Nbd;
        std::size_t B22M = PbM + Mbd;

        std::size_t C11N = PcN;
        std::size_t C11M = PcM;
        std::size_t C12N = PcN;
        std::size_t C12M = PcM + Mbd;
        std::size_t C21N = PcN + Nad;
        std::size_t C21M = PcM;
        std::size_t C22N = PcN + Nad;
        std::size_t C22M = PcM + Mbd;

        //              B  C  BM  SizeofAXX           SizeofBXX           PosOfAXX    PosOfBXX    PosOfCXX
        parlay::par_do([&](){
            parlay::par_do([&](){
                dotDAC(B, C, BM, Nad     , Mad     , Nbd     , Mbd     , A11N, A11M, B11N, B11M, C11N, C11M); // indexing A11, B11, C11
            },[&](){
                dotDAC(B, C, BM, Nad     , Mad     , Nbd     , Mb - Mbd, A11N, A11M, B12N, B12M, C12N, C12M); // indexing A11, B12, C12
            });
        },[&](){
            parlay::par_do([&](){
                dotDAC(B, C, BM, Na - Nad, Mad     , Nbd,      Mbd     , A21N, A21M, B11N, B11M, C21N, C21M); // indexing A21, B11, C21
            },[&](){
                dotDAC(B, C, BM, Na - Nad, Mad     , Nbd, Mb - Mbd     , A21N, A21M, B12N, B12M, C22N, C22M); // indexing A21, B12, C22
            });
        });

        parlay::par_do([&](){
            parlay::par_do([&](){
                dotDAC(B, C, BM, Nad     , Ma - Mad, Nb - Nbd, Mbd     , A12N, A12M, B21N, B21M, C11N, C11M); // indexing A12, B21, C11
            },[&](){
                dotDAC(B, C, BM, Nad     , Ma - Mad, Nb - Nbd, Mb - Mbd, A12N, A12M, B22N, B22M, C12N, C12M); // indexing A12, B22, C12
            });
        },[&](){
            parlay::par_do([&](){
                dotDAC(B, C, BM, Na - Nad, Ma - Mad, Nb - Nbd, Mbd     , A22N, A22M, B21N, B21M, C21N, C21M); // indexing A22, B21, C21
            },[&](){
                dotDAC(B, C, BM, Na - Nad, Ma - Mad, Nb - Nbd, Mb - Mbd, A22N, A22M, B22N, B22M, C22N, C22M); // indexing A22, B22, C22
            });
        });
    }
}

/* Basic inner product */
template<typename P>
void Matrix<P>::dotSeq(const P* B, P* C, std::size_t BM, 
                    std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, 
                    std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM) const {
    for(std::size_t i = 0; i < Na; i++){
        for(std::size_t k = 0; k < Ma; k++){
            for(std::size_t j = 0; j < Mb; j++){
                C[(i+PcN) * BM + (j+PcM)] += data[(i+PaN) * M + (k+PaM)] * B[(k+PbN) * BM + (j+PbM)];
            }
        }
    }
}

/* Returns the transpose of the matrix */
template<typename P>
Matrix<P> Matrix<P>::T() const {
    P* newData = new P[M * N];
    TDAC(newData, N, M, 0, 0);
    return Matrix<P>(M,N,newData);
}

/* Divide and Conquer transpose, allows for parallel execution and better IO efficiency than the basic transpose operation*/
template<typename P>
void Matrix<P>::TDAC(P* D, std::size_t N, std::size_t M, std::size_t PN, std::size_t PM) const {
    if(std::min(N, M) <= 1) {
        TSeq(D, N, M, PN, PM);
    }
    else {
        // size of submatrices
        std::size_t Nd = N / 2;
        std::size_t Md = M / 2;

        // positions of submatrices
        std::size_t D11N = PN;
        std::size_t D11M = PM;
        std::size_t D12N = PN;
        std::size_t D12M = PM + Md;
        std::size_t D21N = PN + Nd;
        std::size_t D21M = PM;
        std::size_t D22N = PN + Nd;
        std::size_t D22M = PM + Md;

        parlay::par_do([&](){
            parlay::par_do([&](){
                TDAC(D, Nd    , Md    , D11N, D11M); // indexing D11
            },[&](){
                TDAC(D, N - Nd, M - Md, D22N, D22M); // indexing D22
            });
        },[&](){
            parlay::par_do([&](){
                TDAC(D, Nd    , M - Md, D12N, D12M); // indexing D12
            },[&](){
                TDAC(D, N - Nd, Md    , D21N, D21M); // indexing D21
            });
        });
    }
}

/* Basic transpose, the data at position i,j is placed into j,i in the new matrix */
template<typename P>
void Matrix<P>::TSeq(P* D, std::size_t N, std::size_t M, std::size_t PN, std::size_t PM) const {
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < M; j++){
            D[(j+PM) * this->N + (i+PN)] = data[(i+PN) * this->M + (j+PM)];
        }
    }
}