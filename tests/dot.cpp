#include <iostream>
#include <vector>
#include <random>
#include "../Matrix.h"

std::vector<std::vector<float>> dot(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(B[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t k = 0; k < A[0].size(); k++){
            for(std::size_t j = 0; j < B[0].size(); j++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

void printMatrix(Matrix<float>& A){
    for(int i = 0; i < A.shape().first; i++){
        for(int j = 0; j < A.shape().second; j++){
            std::cout << A.at(i,j) << " ";
        }
        std::cout << std::endl;
    }
}

void printMatrix(std::vector<std::vector<float>>& A){
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char** argv){
    long int tests = 10;
    char* end;
    if(argc >= 2){
        tests = std::strtol(argv[1], &end, 10);
    }
    if(tests <= 0){
        std::cerr << "Number of tests must be a positive integer" << std::endl;
        return 1;
    }
    auto t = time(0);
    if(argc >= 3){
        t = std::strtol(argv[2], &end, 10);
        if(t <= 0){
            std::cerr << "Seed must be a positive integer" << std::endl;
            return 1;
        }
    }

    std::cout << "Testing inner product " << tests << " times with seed " << t << std::endl;
    srand(t);
    for(long int i = 0; i < tests; i++){
        int M = rand() % 3000 + 1;
        int N = rand() % 3000 + 1;
        int K = rand() % 3000 + 1;

        Matrix<float> A(M,N);
        Matrix<float> B(N,K);


        std::vector<std::vector<float>> Av(M, std::vector<float>(N));
        std::vector<std::vector<float>> Bv(N, std::vector<float>(K));

        for(int a = 0; a < M; a++){
            for(int b = 0; b < N; b++){
                Av[a][b] = rand() % 1000 - 500;
                A.at(a,b) = Av[a][b];
            }
        }
        for(int a = 0; a < N; a++){
            for(int b = 0; b < K; b++){
                Bv[a][b] = rand() % 1000 - 500;
                B.at(a,b) = Bv[a][b];
            }
        }

        Matrix<float> C = A.dot(B);
        std::vector<std::vector<float>> Cv = dot(Av, Bv);

        for(int a = 0; a < M; a++){
            for(int b = 0; b < K; b++){
                if(C.at(a,b) != Cv[a][b]){
                    std::cerr << "Test " << i << " failed at position (" << a << "," << b << ")" << std::endl;
                    std::cout << "Matrix A" << std::endl;
                    printMatrix(A);
                    std::cout << std::endl << "Matrix B" << std::endl;
                    printMatrix(B);
                    std::cout << std::endl << "Matrix C" << std::endl;
                    printMatrix(C);
                    std::cout << std::endl << "Matrix Cv" << std::endl;
                    printMatrix(Cv);
                    return 1;
                }
            }
        }
        std::cout << "Test " << i << " passed " << std::endl;
    }
}