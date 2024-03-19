#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <thread>
#include "parallel.h"

const std::size_t N = 3200;

const int PAR_GRAN = 50;
const int DAC_GRAN = 1000;

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

void DACMult(const int* A, const int* B, int* C, std::size_t Na, std::size_t Ma, std::size_t Nb, std::size_t Mb){
    if(std::min(Na, std::min(Ma, std::min(Nb, Mb))) <= DAC_GRAN) parallelMult(A,B,C,Na,Ma,Nb,Mb);
    else {
        std::size_t Nad = Na / 2;
        std::size_t Mad = Ma / 2;
        std::size_t Nbd = Nb / 2;
        std::size_t Mbd = Mb / 2;
        
        int* B11 = new int[Nbd * Mbd];
        int* B12 = new int[Nbd * (Mb - Mbd)];
        int* B21 = new int[(Nb - Nbd) * Mbd];
        int* B22 = new int[(Nb - Nbd) * (Mb - Mbd)];

        for(std::size_t i = 0; i < Nbd; i++){
            for(std::size_t j = 0; j < Mbd; j++){
                B11[i * Mbd + j] = B[i * Mb + j];
            }
        }
        for(std::size_t i = 0; i < Nbd; i++){
            for(std::size_t j = 0; j < Mb - Mbd; j++){
                B12[i * (Mb - Mbd) + j] = B[i * Mb + (j + Mbd)];
            }
        }
        for(std::size_t i = 0; i < Nb - Nbd; i++){
            for(std::size_t j = 0; j < Mbd; j++){
                B21[i * Mbd + j] = B[((i + Nbd) * Mb) + j];
            }
        }
        for(std::size_t i = 0; i < Nb - Nbd; i++){
            for(std::size_t j = 0; j < Mb - Mbd; j++){
                B22[i * (Mb - Mbd) + j] = B[(i + Nbd) * Mb + (j + Mbd)];
            }
        }

        int* C11 = new int[Nad * Mbd];
        int* C12 = new int[Nad * (Mb - Mbd)];
        int* C21 = new int[(Na - Nad) * Mbd];
        int* C22 = new int[(Na - Nad) * (Mb - Mbd)];

        par_do([&](){
            int* A11 = new int[Nad * Mad];
            for(std::size_t i = 0; i < Nad; i++){
                for(std::size_t j = 0; j < Mad; j++){
                    A11[i * Mad + j] = A[i * Ma + j];
                }
            }
            int* A12 = new int[Nad * (Ma - Mad)];
            for(std::size_t i = 0; i < Nad; i++){
                for(std::size_t j = 0; j < Ma - Mad; j++){
                    A12[i * (Ma - Mad) + j] = A[i * Ma + (j + Mad)];
                }
            }
            par_do([&](){
                int* t1 = new int[Nad * Mbd];
                int* t2 = new int[Nad * Mbd];
                par_do([&](){
                    DACMult(A11, B11, t1, Nad, Mad, Nbd, Mbd); 
                }, [&](){ 
                    DACMult(A12, B21, t2, Nad, Ma - Mad, Nb - Nbd, Mbd); 
                });
                add(t1, t2, C11, Nad, Mbd);
            },[&](){
                int* t1 = new int[Nad * (Mb - Mbd)];
                int* t2 = new int[Nad * (Mb - Mbd)];
                par_do([&](){ 
                    DACMult(A11, B12, t1, Nad, Mad, Nbd, Mb - Mbd); 
                }, [&](){ 
                    DACMult(A12, B22, t2, Nad, Ma - Mad, Nb - Nbd, Mb - Mbd); 
                });
                add(t1, t2, C12, Nad, Mb - Mbd);
            });
        },[&](){
            int* A21 = new int[(Na - Nad) * Mad];
            for(std::size_t i = 0; i < Na - Nad; i++){
                for(std::size_t j = 0; j < Mad; j++){
                    A21[i * Mad + j] = A[(i + Nad) * Ma + j];
                }
            }
            int* A22 = new int[(Na - Nad) * (Ma - Mad)];
            for(std::size_t i = 0; i < Na - Nad; i++){
                for(std::size_t j = 0; j < Ma - Mad; j++){
                    A22[i * (Ma - Mad) + j] = A[(i + Nad) * Ma + (j + Mad)];
                }
            }
            par_do([&](){
                int* t1 = new int[(Na - Nad) * Mbd];
                int* t2 = new int[(Na - Nad) * Mbd];
                par_do([&](){ 
                    DACMult(A21, B11, t1, Na - Nad, Mad, Nbd, Mbd); 
                }, [&](){ 
                    DACMult(A22, B21, t2, Na - Nad, Ma - Mad, Nb - Nbd, Mbd); 
                });
                add(t1, t2, C21, Na - Nad, Mbd);
            },[&](){
                int* t1 = new int[(Na - Nad) * (Mb - Mbd)];
                int* t2 = new int[(Na - Nad) * (Mb - Mbd)];
                par_do([&](){ 
                    DACMult(A21, B12, t1, Na - Nad, Mad, Nbd, Mb - Mbd); 
                }, [&](){ 
                    DACMult(A22, B22, t2, Na - Nad, Ma - Mad, Nb - Nbd, Mb - Mbd); 
                });
                add(t1, t2, C22, Na - Nad, Mb - Mbd);
            });
        });

        for(std::size_t i = 0; i < Nad; i++){
            for(std::size_t j = 0; j < Mbd; j++){
                C[i * Mb + j] = C11[i * Mbd + j];
            }
        }
        for(std::size_t i = 0; i < Nad; i++){
            for(std::size_t j = 0; j < Mb-Mbd; j++){
                C[i * Mb + (j + Mbd)] = C12[i * (Mb-Mbd) + j];
            }
        }
        for(std::size_t i = 0; i < Na - Nad; i++){
            for(std::size_t j = 0; j < Mbd; j++){
                C[(i + Nad) * Mb + j] = C21[i * Mbd + j];
            }
        }
        for(std::size_t i = 0; i < Na - Nad; i++){
            for(std::size_t j = 0; j < Mb - Mbd; j++){
                C[(i + Nad) * Mb + (j + Mbd)] = C22[i * (Mb - Mbd) + j];
            }
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
        }
    }

    std::vector<double> times(3);

    for(int i = 0; i < 4; i++){
        std::cout << "Round " << i + 1 << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        sequentialMult(A, B, C, N, N, N, N);
        times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        std::cout << "finished sequential" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        parallelMult(A, B, C, N, N, N, N);
        times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        std::cout << "finished parallel" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        DACMult(A, B, C, N, N, N, N);
        times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        std::cout << "finished DAC" << std::endl;
    }
    std::cout << "Sequential: " << times[0] / 4 << "ms" << std::endl;
    std::cout << "Parallel: " << times[1] / 4 << "ms" << std::endl;
    std::cout << "Divide and Conquer: " << times[2] / 4 << "ms" << std::endl;
}