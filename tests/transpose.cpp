#include <iostream>
#include <vector>
#include <random>
#include "../Matrix.h"

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>>& A){
    std::vector<std::vector<double>> C(A[0].size(), std::vector<double>(A.size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[j][i] = A[i][j];
        }
    }
    return C;
}

void printMatrix(Matrix<double>& A){
    for(int i = 0; i < A.shape().first; i++){
        for(int j = 0; j < A.shape().second; j++){
            std::cout << A.at(i,j) << " ";
        }
        std::cout << std::endl;
    }
}

void printMatrix(std::vector<std::vector<double>>& A){
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

    std::cout << "Testing Transpose " << tests << " times with seed " << t << std::endl;
    srand(t);
    for(long int i = 0; i < tests; i++){
        int N = rand() % 3000 + 1;
        int M = rand() % 3000 + 1;

        Matrix<double> A(N,M);

        std::vector<std::vector<double>> Av(N, std::vector<double>(M));

        for(int a = 0; a < N; a++){
            for(int b = 0; b < M; b++){
                Av[a][b] = rand() % 10;
                A.at(a,b) = Av[a][b];
            }
        }

        Matrix<double> B = A.T();
        std::vector<std::vector<double>> Bv = transpose(Av);

        for(int a = 0; a < M; a++){
            for(int b = 0; b < N; b++){
                if(B.at(a,b) != Bv[a][b]){
                    std::cerr << "Test " << i << " failed at position (" << a << "," << b << ")" << std::endl;
                    std::cout << "Matrix" << std::endl;
                    printMatrix(A);
                    std::cout << "Transpose" << std::endl;
                    printMatrix(B);
                    std::cout << "TransposeV" << std::endl;
                    printMatrix(Bv);
                    return 1;
                }
            }
        }
        std::cout << "Test " << i << " passed " << std::endl;
    }
}