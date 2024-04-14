#pragma once
#include <cstddef>
#include <utility>
#include <string>
#include <cmath>
#include <limits>
#include "parallel.h"

const std::size_t GRANULARITY = 30000;
// const std::size_t GRANULARITY = 300;

template<typename P>
class Matrix{
    public:
    Matrix() : M(0), N(0) {} // default constructor
    Matrix(std::size_t,std::size_t);
    Matrix(std::size_t,std::size_t,P);
    Matrix(std::size_t m, std::size_t n, std::vector<P>& d);
    ~Matrix();
    std::pair<std::size_t,std::size_t> shape() const;
    void reshape(std::size_t,std::size_t);
    Matrix<P> rRange(std::size_t, std::size_t) const;
    Matrix<P> cRange(std::size_t, std::size_t) const;
    Matrix<P> dot(const Matrix<P>&) const;
    Matrix<P> add(const Matrix<P>&) const;
    Matrix<P> sub(const Matrix<P>&) const;
    Matrix<P> mult(const Matrix<P>&) const;
    Matrix<P> div(const Matrix<P>&) const;
    Matrix<P> log() const;
    Matrix<P> exp() const;
    Matrix<P> tanh() const;
    std::pair<std::size_t,std::size_t> argMax() const;
    Matrix<P> argMax(bool) const;
    P max() const;
    Matrix<P> max(bool) const;
    P sum() const;
    Matrix<P> sum(bool) const;
    Matrix<P> T() const;
    P& at(std::size_t,std::size_t);
    P at(std::size_t, std::size_t) const;
    template<typename U> Matrix<P> add(const U&) const;
    template<typename U> Matrix<P> sub(const U&) const;
    template<typename U> Matrix<P> subL(const U&) const;
    template<typename U> Matrix<P> mult(const U&) const;
    template<typename U> Matrix<P> div(const U&) const;
    template<typename U> Matrix<P> divL(const U&) const;
    std::string strShape() const;
    bool operator==(const Matrix<P>&);
    private:
    std::size_t M;
    std::size_t N;
    std::vector<P> data;
    void dotDAC(const std::vector<P>&, std::vector<P>&, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t) const;
    void dotSeq(const std::vector<P>&, std::vector<P>&, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t) const;
    void TDAC(std::vector<P>&, std::size_t, std::size_t, std::size_t, std::size_t) const;
    void TSeq(std::vector<P>&, std::size_t, std::size_t, std::size_t, std::size_t) const;
};

/* Constructor */
template<typename P>
Matrix<P>::Matrix(std::size_t m, std::size_t n, P v) : M(m), N(n){
    if(!M || !N) throw std::runtime_error("Invalid matrix dimensions: " + strShape());
    data = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M * N; i++){
        data[i] = v;
    }
}

/* Constructor with data */
template<typename P>
Matrix<P>::Matrix(std::size_t m, std::size_t n, std::vector<P>& d) : M(m), N(n), data(d){
    if(M*N != data.size())
        throw std::runtime_error("Cannot Construct matrix with shape " + strShape + " and vector of size " data.size());
}

/* Destructor */
template<typename P>
Matrix<P>::~Matrix(){}

/* Constructor with shape */
template<typename P>
Matrix<P>::Matrix(std::size_t m, std::size_t n) : M(m), N(n){
    if(M <= 0 || N <= 0) throw std::runtime_error("Invalid matrix dimensions: " + strShape());
    data = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M * N; i++){
        data[i] = P();
    }
}

/* Returns a pair with the number of rows and the number of columns in the matrix */
template<typename P>
std::pair<std::size_t,std::size_t> Matrix<P>::shape() const {
    return std::make_pair(M,N);
}

/* Update M and N if new values are compatible */
template<typename P>
void Matrix<P>::reshape(std::size_t m, std::size_t n){
    if(M*N != m*n) throw std::runtime_error("Incompatible matrix dimensions for reshape: " + strShape() + " (" + std::to_string(m) + "," + std::to_string(n) + ")");
    M = m;
    N = n;
}

/* Returns a string representation of the shape of the matrix: "(M,N)" */
template<typename P>
std::string Matrix<P>::strShape() const {
    return "(" + std::to_string(M) + "," + std::to_string(N) + ")";
}

/* Returns a reference to the value in the matrix at the given coordinates */
template<typename P>
P& Matrix<P>::at(std::size_t i, std::size_t j) {
    return data[i * N + j];
}

/* Returns the value in the matrix at the given coordinates */
template<typename P>
P Matrix<P>::at(std::size_t i, std::size_t j) const {
    return data[i * N + j];
}

/* Returns the submatrix of rows in the given range(non inclusive)*/
template<typename P>
Matrix<P> Matrix<P>::rRange(std::size_t from, std::size_t to) const {
    if(to <= from)
        return Matrix<P>();
    std::vector<P> nData((to-from) * N);
    for(std::size_t i = 0; i < to-from; i++){
        for(std::size_t j = 0; j < N; j++){
            nData[i * N + j] = data[(i+from) * N + j];
        }
    }
    return Matrix<P>(to-from, N, nData);
}

/* Returns the submatrix of columns in the given range(non inclusive)*/
template<typename P>
Matrix<P> Matrix<P>::cRange(std::size_t from, std::size_t to) const {
    if(to <= from)
        return Matrix<P>();
    std::vector<P> nData(M * (to-from));
    for(std::size_t i = 0; i < M; i++){
        for(std::size_t j = 0; j < (to-from); j++){
            nData[i * (to-from) + j] = data[i * N + j+from];
        }
    }
    return Matrix<P>(M, to-from, nData);
}

/* Computes the dot product of two matrices, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::dot(const Matrix<P>& B) const {
    if(N != B.M) throw std::runtime_error("Incompatible matrix dimensions for dot product: " + strShape() + " " + B.strShape());
    Matrix C(M, B.N);
    dotDAC(B.data, C.data, B.N, M, N, B.M, B.N, 0, 0, 0, 0, 0, 0);
    return C;
}

/*
Element-wise addition of two matrices, returns a matrix 
If the matrices have different dimensions, the smaller matrix is broadcasted to the larger matrix
*/
template<typename P>
Matrix<P> Matrix<P>::add(const Matrix<P>& B) const {
    if(M != B.M && N != B.N) throw std::runtime_error("Incompatible matrix dimensions for addition: " + strShape() + " " + B.strShape());
    if(M == B.M && N == B.N){
        std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
        for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] + B.data[i];
        return Matrix<P>(M, N, Cdata);
    }
    std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
    for(std::size_t i = 0; i < std::max(M,B.M); i++){
        for(std::size_t j = 0; j < std::max(N,B.N); j++){
            Cdata[i * std::max(N,B.N) + j] = data[(i%M) * N + (j%N)] + B.data[(i%B.M) * B.N + (j%B.N)];
        }
    }
    return Matrix<P>(std::max(M,B.M),std::max(N,B.N),Cdata);
}
/* Element-wise addition of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::add(const U& v) const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] + v;
    return Matrix<P>(M,N,Cdata);
}
/*
Element-wise multiplication of two matrices, returns a matrix 
If the matrices have different dimensions, the smaller matrix is broadcasted to the larger matrix
*/
template<typename P>
Matrix<P> Matrix<P>::mult(const Matrix<P>& B) const {
    if(M != B.M && N != B.N) throw std::runtime_error("Incompatible matrix dimensions for multiplication: " + strShape() + " " + B.strShape());
    if(M == B.M && N == B.N){
        std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
        for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] * B.data[i];
        return Matrix<P>(M, N, Cdata);
    }
    std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
    for(std::size_t i = 0; i < std::max(M,B.M); i++){
        for(std::size_t j = 0; j < std::max(N,B.N); j++){
            Cdata[i * std::max(N,B.N) + j] = data[(i%M) * N + (j%N)] * B.data[(i%B.M) * B.N + (j%B.N)];
        }
    }
    return Matrix<P>(std::max(M,B.M),std::max(N,B.N),Cdata);
}
/* Element-wise multiplication of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::mult(const U& v) const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] * v;
    return Matrix<P>(M,N,Cdata);
}
/*
Element-wise subtraction of two matrices, returns a matrix 
If the matrices have different dimensions, the smaller matrix is broadcasted to the larger matrix
*/
template<typename P>
Matrix<P> Matrix<P>::sub(const Matrix<P>& B) const {
    if(M != B.M && N != B.N) throw std::runtime_error("Incompatible matrix dimensions for subtraction: " + strShape() + " " + B.strShape());
    if(M == B.M && N == B.N){
        std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
        for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] - B.data[i];
        return Matrix<P>(M, N, Cdata);
    }
    std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
    for(std::size_t i = 0; i < std::max(M,B.M); i++){
        for(std::size_t j = 0; j < std::max(N,B.N); j++){
            Cdata[i * std::max(N,B.N) + j] = data[(i%M) * N + (j%N)] - B.data[(i%B.M) * B.N + (j%B.N)];
        }
    }
    return Matrix<P>(std::max(M,B.M),std::max(N,B.N),Cdata);
}
/* Element-wise subtraction of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::sub(const U& v) const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] - v;
    return Matrix<P>(M,N,Cdata);
}
/* Element-wise subtraction of a scalar on the left, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::subL(const U& v) const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = v - data[i];
    return Matrix<P>(M,N,Cdata);
}
/*
Element-wise division of two matrices, returns a matrix 
If the matrices have different dimensions, the smaller matrix is broadcasted to the larger matrix
*/
template<typename P>
Matrix<P> Matrix<P>::div(const Matrix<P>& B) const {
    if(M != B.M && N != B.N) throw std::runtime_error("Incompatible matrix dimensions for division: " + strShape() + " " + B.strShape());
    if(M == B.M && N == B.N){
        std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
        for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] / B.data[i];
        return Matrix<P>(M, N, Cdata);
    }
    std::vector<P> Cdata = std::vector<P>(std::max(M,B.M)*std::max(N,B.N));
    for(std::size_t i = 0; i < std::max(M,B.M); i++){
        for(std::size_t j = 0; j < std::max(N,B.N); j++){
            Cdata[i * std::max(N,B.N) + j] = data[(i%M) * N + (j%N)] / B.data[(i%B.M) * B.N + (j%B.N)];
        }
    }
    return Matrix<P>(std::max(M,B.M),std::max(N,B.N),Cdata);
}
/* Element-wise division of a scalar, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::div(const U& v) const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = data[i] / v;
    return Matrix<P>(M,N,Cdata);
}
/* Element-wise division of a scalar on the left, returns a matrix */
template<typename P>
template<typename U>
Matrix<P> Matrix<P>::divL(const U& v) const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = v / data[i];
    return Matrix<P>(M,N,Cdata);
}
/* Element-wise logarithm of the matrix, returns a matrix*/
template<typename P>
Matrix<P> Matrix<P>::log() const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = std::log(data[i]);
    return Matrix<P>(M,N,Cdata);
}
/* Sum of all elements in the matrix */
template<typename P>
P Matrix<P>::sum() const {
    P s = P();
    for(std::size_t i = 0; i < M*N; i++) s += data[i];
    return s;
}
/* Returns a matrix of the sum of each row if 0 and each column if 1 */
template<typename P>
Matrix<P> Matrix<P>::sum(bool r) const {
    std::size_t m = 0;
    std::size_t n = 0;
    if(!r)
        m = M;
    else
        n = N;
    std::vector<P> Cdata(m+n, 0);
    for(std::size_t i = 0; i < M; i++){
        for(std::size_t j = 0; j < N; j++){
            Cdata[i*(!r)+j*r] += data[i*N + j];
        }
    }
    if(m) n = 1;
    else m = 1;
    return Matrix<P>(m,n,Cdata);
}
/* Element-wise exponentiation of the matrix, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::exp() const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = std::exp(data[i]);
    return Matrix<P>(M,N,Cdata);
}

/* Element-wise tanh of the matrix, returns a matrix */
template<typename P>
Matrix<P> Matrix<P>::tanh() const {
    std::vector<P> Cdata = std::vector<P>(M*N);
    for(std::size_t i = 0; i < M*N; i++) Cdata[i] = std::tanh(data[i]);
    return Matrix<P>(M,N,Cdata);
}

/* Returns the position of the maximum element in the matrix */
template<typename P>
std::pair<std::size_t,std::size_t> Matrix<P>::argMax() const {
    std::size_t p1 = 0;
    std::size_t p2 = 0;
    for(std::size_t i = 0; i < M; i++){
        for(std::size_t j = 0; j < N; j++){
            if(data[p1*N+p2] < data[i*N+j]){
                p1=i;
                p2=j;
            }
        }
    }
    return std::make_pair(p1,p2);
}

/* Returns the position of the maximum element in each row if 0 and each column if 1 */
template<typename P>
Matrix<P> Matrix<P>::argMax(bool r) const {
    std::size_t m = 0;
    std::size_t n = 0;
    if(!r)
        m = M;
    else
        n = N;
    std::vector<P> Cdata(m+n,0);
    for(std::size_t i = 0; i < M; i++){
        for(std::size_t j = 0; j < N; j++){
            int CdataPosition = i*(!r)+j*r; // position in C constrained to an axis
            int dataPositionWRTCdata = (i*(!r)+Cdata[CdataPosition]*r) * N + (j*r+Cdata[CdataPosition]*(!r)); // position in data given the data in Cdata constrained to an axis
            int dataPosition = i*N + j; // current position in the data
            if(data[dataPositionWRTCdata] < data[dataPosition])
                Cdata[CdataPosition] = i*r+j*(!r);
        }
    }
    if(m) n = 1;
    else m = 1;
    return Matrix<P>(m,n,Cdata);
}

/* Returns the maximum element in the matrix */
template<typename P>
P Matrix<P>::max() const {
    P m = std::numeric_limits<P>::lowest();
    for(std::size_t i = 0; i < M*N; i++) m = std::max(m, data[i]);
    return m;
}

/* Returns a matrix of the maximum element in each row if 0 and each column if 1 */
template<typename P>
Matrix<P> Matrix<P>::max(bool r) const{
    Matrix<P> vals = argMax(r);
    for(std::size_t p = 0; p < M*(!r) + N*r; p++)
        vals.data[p] = data[(vals.data[p]*r+p*(!r))*N+(vals.data[p]*(!r)+p*r)];
    return vals;
}

/* Divide and Conquer dot product, allows for parallel execution and better IO efficiency than the basic inner product operation */
template<typename P>
void Matrix<P>::dotDAC(const std::vector<P>& B, std::vector<P>& C, std::size_t BN, 
             std::size_t Ma, std::size_t Na, std::size_t Mb, std::size_t Nb, 
             std::size_t PaM, std::size_t PaN, std::size_t PbM, std::size_t PbN, std::size_t PcM, std::size_t PcN) const {
    if(Ma*Na <= GRANULARITY || Mb*Nb <= GRANULARITY) {
    // if(std::min(std::min(Ma, Na), std::min(Mb, Nb)) <= GRANULARITY){
        dotSeq(B,C,BN,Ma,Na,Mb,Nb,PaM,PaN,PbM,PbN,PcM,PcN);
    }
    else {
        // size of submatrices
        std::size_t Mad = Ma / 2;
        std::size_t Nad = Na / 2;
        std::size_t Mbd = Mb / 2;
        std::size_t Nbd = Nb / 2;

        // positions of submatrices
        std::size_t A11M = PaM;
        std::size_t A11N = PaN;
        std::size_t A12M = PaM;
        std::size_t A12N = PaN + Nad;
        std::size_t A21M = PaM + Mad;
        std::size_t A21N = PaN;
        std::size_t A22M = PaM + Mad;
        std::size_t A22N = PaN + Nad;

        std::size_t B11M = PbM;
        std::size_t B11N = PbN;
        std::size_t B12M = PbM;
        std::size_t B12N = PbN + Nbd;
        std::size_t B21M = PbM + Mbd;
        std::size_t B21N = PbN;
        std::size_t B22M = PbM + Mbd;
        std::size_t B22N = PbN + Nbd;

        std::size_t C11M = PcM;
        std::size_t C11N = PcN;
        std::size_t C12M = PcM;
        std::size_t C12N = PcN + Nbd;
        std::size_t C21M = PcM + Mad;
        std::size_t C21N = PcN;
        std::size_t C22M = PcM + Mad;
        std::size_t C22N = PcN + Nbd;

        //             B  C  BN  SizeofAXX           SizeofBXX           PosOfAXX    PosOfBXX    PosOfCXX
        parlay::par_do([&](){
            parlay::par_do([&](){
                dotDAC(B, C, BN, Mad     , Nad     , Mbd     , Nbd     , A11M, A11N, B11M, B11N, C11M, C11N); // indexing A11, B11, C11
            },[&](){
                dotDAC(B, C, BN, Mad     , Nad     , Mbd     , Nb - Nbd, A11M, A11N, B12M, B12N, C12M, C12N); // indexing A11, B12, C12
            });
        },[&](){
            parlay::par_do([&](){
                dotDAC(B, C, BN, Ma - Mad, Nad     , Mbd     ,      Nbd, A21M, A21N, B11M, B11N, C21M, C21N); // indexing A21, B11, C21
            },[&](){
                dotDAC(B, C, BN, Ma - Mad, Nad     , Mbd     , Nb - Nbd, A21M, A21N, B12M, B12N, C22M, C22N); // indexing A21, B12, C22
            });
        });

        parlay::par_do([&](){
            parlay::par_do([&](){
                dotDAC(B, C, BN, Mad     , Na - Nad, Mb - Mbd, Nbd     , A12M, A12N, B21M, B21N, C11M, C11N); // indexing A12, B21, C11
            },[&](){
                dotDAC(B, C, BN, Mad     , Na - Nad, Mb - Mbd, Nb - Nbd, A12M, A12N, B22M, B22N, C12M, C12N); // indexing A12, B22, C12
            });
        },[&](){
            parlay::par_do([&](){
                dotDAC(B, C, BN, Ma - Mad, Na - Nad, Mb - Mbd, Nbd     , A22M, A22N, B21M, B21N, C21M, C21N); // indexing A22, B21, C21
            },[&](){
                dotDAC(B, C, BN, Ma - Mad, Na - Nad, Mb - Mbd, Nb - Nbd, A22M, A22N, B22M, B22N, C22M, C22N); // indexing A22, B22, C22
            });
        });
    }
}

/* Basic dot product */
template<typename P>
void Matrix<P>::dotSeq(const std::vector<P>& B, std::vector<P>& C, std::size_t BN, 
                    std::size_t Ma, std::size_t Na, std::size_t Mb, std::size_t Nb, 
                    std::size_t PaM, std::size_t PaN, std::size_t PbM, std::size_t PbN, std::size_t PcM, std::size_t PcN) const {
    for(std::size_t i = 0; i < Ma; i++){
        for(std::size_t k = 0; k < Na; k++){
            for(std::size_t j = 0; j < Nb; j++){
                C[(i+PcM) * BN + (j+PcN)] += data[(i+PaM) * N + (k+PaN)] * B[(k+PbM) * BN + (j+PbN)];
            }
        }
    }
}

/* Returns the transpose of the matrix */
template<typename P>
Matrix<P> Matrix<P>::T() const {
    std::vector<P> newData = std::vector<P>(M*N);
    TDAC(newData, M, N, 0, 0);
    return Matrix<P>(N,M,newData);
}

/* Divide and Conquer transpose, allows for parallel execution and better IO efficiency than the basic transpose operation */
template<typename P>
void Matrix<P>::TDAC(std::vector<P>& D, std::size_t M, std::size_t N, std::size_t PM, std::size_t PN) const {
    if(M*N <= GRANULARITY) {
    // if(std::min(M,N) <= GRANULARITY){
        TSeq(D, M, N, PM, PN);
    }
    else {
        // size of submatrices
        std::size_t Md = M / 2;
        std::size_t Nd = N / 2;

        // positions of submatrices
        std::size_t D11M = PM;
        std::size_t D11N = PN;
        std::size_t D12M = PM;
        std::size_t D12N = PN + Nd;
        std::size_t D21M = PM + Md;
        std::size_t D21N = PN;
        std::size_t D22M = PM + Md;
        std::size_t D22N = PN + Nd;

        parlay::par_do([&](){
            parlay::par_do([&](){
                TDAC(D, Md    , Nd    , D11M, D11N); // indexing D11
            },[&](){
                TDAC(D, M - Md, N - Nd, D22M, D22N); // indexing D22
            });
        },[&](){
            parlay::par_do([&](){
                TDAC(D, Md    , N - Nd, D12M, D12N); // indexing D12
            },[&](){
                TDAC(D, M - Md, Nd    , D21M, D21N); // indexing D21
            });
        });
    }
}

/* Basic transpose, the data at position i,j is placed into j,i in the new matrix */
template<typename P>
void Matrix<P>::TSeq(std::vector<P>& D, std::size_t M, std::size_t N, std::size_t PM, std::size_t PN) const {
    for(std::size_t i = 0; i < M; i++){
        for(std::size_t j = 0; j < N; j++){
            D[(j+PN) * this->M + (i+PM)] = data[(i+PM) * this->N + (j+PN)];
        }
    }
}

/* Equality operator */
template<typename P>
bool Matrix<P>::operator==(const Matrix<P>& rhs){
    if(M != rhs.M || N != rhs.N) return false;
    if(data != rhs.data) return false;
    return true;
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
    return rhs.subL(lhs);
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
/* Element-wise division of matrices */
template<typename P>
inline Matrix<P> operator/(const Matrix<P>& lhs, const Matrix<P>& rhs){
    return lhs.div(rhs);
}
/* Element-wise division of a scalar on the left hand side */
template<typename P, typename U>
inline Matrix<P> operator/(const U& lhs, const Matrix<P>& rhs){
    return rhs.divL(lhs);
}
/* Element-wise division of a scalar on the right hand side */
template<typename P, typename U>
inline Matrix<P> operator/(const Matrix<P>& lhs, const U& rhs){
    return lhs.div(rhs);
}