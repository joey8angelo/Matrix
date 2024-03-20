#include "Matrix.h"
#include "parallel.h"

Matrix::Matrix(int n, int m) : N(n), M(m){
    data = new double[N * M];
    for(std::size_t i = 0; i < N * M; i++){
        data[i] = 0;
    }
}

Matrix::~Matrix(){
    delete[] data;
}

std::pair<int,int> Matrix::shape(){
    return std::make_pair(N,M);
}

void Matrix::reshape(int n, int m){
    return;
}

const Matrix Matrix::mult(const Matrix& B){
    if(M != B.N) throw "Incompatible matrix dimensions" + std::to_string(M) + " != " + std::to_string(B.N);
    Matrix C(N, B.M);
    DACMult(B.data, C.data, B.M, N, M, B.N, B.M, 0, 0, 0, 0, 0, 0);
    return C;
}

double& Matrix::at(int i, int j){
    return data[i * M + j];
}

const void Matrix::DACMult(const double* B, double* C, std::size_t BM, 
             std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, 
             std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM){
    if(std::min(Na, std::min(Ma, std::min(Nb, Mb))) <= DAC_GRANULAR) {
        multPar(B,C,BM,Na,Ma,Nb,Mb,PaN,PaM,PbN,PbM,PcN,PcM);
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

        //                  B  C  BM  SizeofAXX           SizeofBXX           PosOfAXX    PosOfBXX    PosOfCXX
        parlay::par_do([&](){
            parlay::par_do([&](){
                parlay::par_do([&](){
                    DACMult(B, C, BM, Nad     , Mad     , Nbd     , Mbd     , A11N, A11M, B11N, B11M, C11N, C11M); // indexing A11, B11, C11
                },[&](){
                    DACMult(B, C, BM, Nad     , Mad     , Nbd     , Mb - Mbd, A11N, A11M, B12N, B12M, C12N, C12M); // indexing A11, B12, C12
                });
            },[&](){
                parlay::par_do([&](){
                    DACMult(B, C, BM, Na - Nad, Mad     , Nbd,      Mbd     , A21N, A21M, B11N, B11M, C21N, C21M); // indexing A21, B11, C21
                },[&](){
                    DACMult(B, C, BM, Na - Nad, Mad     , Nbd, Mb - Mbd     , A21N, A21M, B12N, B12M, C22N, C22M); // indexing A21, B12, C22
                });
            });
        },[&](){
            parlay::par_do([&](){
                parlay::par_do([&](){
                    DACMult(B, C, BM, Nad     , Ma - Mad, Nb - Nbd, Mbd     , A12N, A12M, B21N, B21M, C11N, C11M); // indexing A12, B21, C11
                },[&](){
                    DACMult(B, C, BM, Nad     , Ma - Mad, Nb - Nbd, Mb - Mbd, A12N, A12M, B22N, B22M, C12N, C12M); // indexing A12, B22, C12
                });
            },[&](){
                parlay::par_do([&](){
                    DACMult(B, C, BM, Na - Nad, Ma - Mad, Nb - Nbd, Mbd     , A22N, A22M, B21N, B21M, C21N, C21M); // indexing A22, B21, C21
                },[&](){
                    DACMult(B, C, BM, Na - Nad, Ma - Mad, Nb - Nbd, Mb - Mbd, A22N, A22M, B22N, B22M, C22N, C22M); // indexing A22, B22, C22
                });
            });
        });
    }
}

const void Matrix::multPar(const double* B, double* C, std::size_t BM, 
                  std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, 
                  std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM){
    if(std::min(Na, std::min(Ma, std::min(Nb, Mb))) <= PAR_GRANULAR) {
        multSeq(B,C,BM,Na,Ma,Nb,Mb,PaN,PaM,PbN,PbM,PcN,PcM);
    }
    else{
        parlay::parallel_for(0, Na, [&](int i){
            for(std::size_t k = 0; k < Ma; k++){
                for(std::size_t j = 0; j < Mb; j++){
                    C[(i+PcN) * BM + (j+PcM)] += data[(i+PaN) * M + (k+PaM)] * B[(k+PbN) * BM + (j+PbM)];
                }
            }
        });
    }
}

const void Matrix::multSeq(const double* B, double* C, std::size_t BM, 
                    std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, 
                    std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM){
    for(std::size_t i = 0; i < Na; i++){
        for(std::size_t k = 0; k < Ma; k++){
            for(std::size_t j = 0; j < Mb; j++){
                C[(i+PcN) * BM + (j+PcM)] += data[(i+PaN) * M + (k+PaM)] * B[(k+PbN) * BM + (j+PbM)];
            }
        }
    }
}