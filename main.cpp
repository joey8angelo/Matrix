#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <thread>
#include "parallel.h"

const std::size_t N = 3200;

const int PAR_GRAN = 50;
const int DAC_GRAN = 50;

using namespace parlay;

void sequentialMult(const int* A, const int* B, int* C, std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb){
    for(std::size_t i = 0; i < Na; i++){
        for(std::size_t k = 0; k < Ma; k++){
            for(std::size_t j = 0; j < Mb; j++){
                C[i * Mb + j] = A[i * Ma + k] * B[k * Mb + j];
            }
        }
    }
}

void parallelMult(const int* A, const int* B, int* C, std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb){
    if(std::min(Na, std::min(Ma, std::min(Nb, Mb))) <= PAR_GRAN) return sequentialMult(A, B, C, Na, Ma, Nb, Mb);
    parallel_for(0, Na, [&](int i){
        for(std::size_t k = 0; k < Ma; k++){
            for(std::size_t j = 0; j < Mb; j++){
                C[i * Mb + j] += A[i * Ma + k] * B[k * Mb + j];
            }
        }
    });
}

void add(const int* A, const int* B, int* C, std::size_t N, std::size_t M){
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < M; j++){
            C[i * M + j] = A[i * M + j] + B[i * M + j];
        }
    }
}

void multAt(const int* A, const int* B, int* C, std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM){
    for(std::size_t i = 0; i < Na; i++){
        for(std::size_t k = 0; k < Ma; k++){
            for(std::size_t j = 0; j < Mb; j++){
                C[(i+PcN) * N + (j+PcM)] += A[(i+PaN) * N + (k+PaM)] * B[(k+PbN) * N + (j+PbM)];
            }
        }
    }
}

void DACMult(const int* A, const int* B, int* C, std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb, std::size_t PaN, std::size_t PaM, std::size_t PbN, std::size_t PbM, std::size_t PcN, std::size_t PcM){
    if(std::min(Na, std::min(Ma, std::min(Nb, Mb))) <= DAC_GRAN) multAt(A,B,C,Na,Ma,Nb,Mb,PaN,PaM,PbN,PbM,PcN,PcM);
    else {
        std::size_t Nad = Na / 2;
        std::size_t Mad = Ma / 2;
        std::size_t Nbd = Nb / 2;
        std::size_t Mbd = Mb / 2;

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

        par_do([&](){
            par_do([&](){
                par_do([&](){
                    DACMult(A, B, C, Nad, Mad     , Nbd     , Mbd, A11N, A11M, B11N, B11M, C11N, C11M); // indexing A11, B11, C11
                },[&](){
                    DACMult(A, B, C, Nad, Mad     ,      Nbd, Mb - Mbd, A11N, A11M, B12N, B12M, C12N, C12M); // indexing A11, B12, C12
                });
            },[&](){
                par_do([&](){
                    DACMult(A, B, C, Na - Nad, Mad     , Nbd     , Mbd, A21N, A21M, B11N, B11M, C21N, C21M); // indexing A21, B11, C21
                },[&](){
                    DACMult(A, B, C, Na - Nad, Mad     , Nbd     , Mb - Mbd, A21N, A21M, B12N, B12M, C22N, C22M); // indexing A21, B12, C22
                });
            });
        },[&](){
            par_do([&](){
                par_do([&](){
                    DACMult(A, B, C, Nad, Ma - Mad, Nb - Nbd, Mbd, A12N, A12M, B21N, B21M, C11N, C11M); // indexing A12, B21, C11
                },[&](){
                    DACMult(A, B, C, Nad, Ma - Mad, Nb - Nbd, Mb - Mbd, A12N, A12M, B22N, B22M, C12N, C12M); // indexing A12, B22, C12
                });
            },[&](){
                par_do([&](){
                    DACMult(A, B, C, Na - Nad, Ma - Mad, Nb - Nbd, Mbd, A22N, A22M, B21N, B21M, C21N, C21M); // indexing A22, B21, C21
                },[&](){
                    DACMult(A, B, C, Na - Nad, Ma - Mad, Nb - Nbd, Mb - Mbd, A22N, A22M, B22N, B22M, C22N, C22M); // indexing A22, B22, C22
                });
            });
        });

    }
}

void zeroOut(int* A, std::size_t N, std::size_t M){
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < M; j++){
            A[i * M + j] = 0;
        }
    }
}

int main(){
    int* A = new int[N*N];
    int* B = new int[N*N];
    int* C = new int[N*N];
    srand(0);
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < N; j++){
            A[i*N + j] = rand() % 10;
            B[i*N + j] = rand() % 10;
            C[i*N + j] = 0;
        }
    }

    std::vector<double> times(3);

    for(int i = 0; i < 4; i++){
        std::cout << "Round " << i + 1 << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        sequentialMult(A, B, C, N, N, N, N);
        times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        std::cout << "finished sequential" << std::endl;
        zeroOut(C, N, N);
        start = std::chrono::high_resolution_clock::now();
        parallelMult(A, B, C, N, N, N, N);
        times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        std::cout << "finished parallel" << std::endl;
        zeroOut(C, N, N);
        start = std::chrono::high_resolution_clock::now();
        DACMult(A, B, C, N, N, N, N, 0, 0, 0, 0, 0, 0);
        times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        std::cout << "finished DAC" << std::endl << std::endl;
        zeroOut(C, N, N);
    }
    std::cout << "Sequential: " << times[0] / 4 << "ms" << std::endl;
    std::cout << "Parallel: " << times[1] / 4 << "ms" << std::endl;
    std::cout << "Divide and Conquer: " << times[2] / 4 << "ms" << std::endl;
}