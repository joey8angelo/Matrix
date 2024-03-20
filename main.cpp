#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include "Matrix.h"

int main(){
    const std::size_t N = 3200;
    Matrix A(N,N);
    Matrix B(N,N);
    srand(0);
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < N; j++){
            A.at(i,j) = rand() % 10;
            B.at(i,j) = rand() % 10;
        }
    }

    // for(std::size_t i = 0; i < N; i++){
    //     for(std::size_t j = 0; j < N; j++){
    //         std::cout << A.at(i,j) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    // for(std::size_t i = 0; i < N; i++){
    //     for(std::size_t j = 0; j < N; j++){
    //         std::cout << B.at(i,j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    Matrix Ct = A.mult(B);
    // int n = Ct.shape().first;
    // int m = Ct.shape().second;
    // for(std::size_t i = 0; i < n; i++){
    //     for(std::size_t j = 0; j < m; j++){
    //         std::cout << Ct.at(i,j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    int time = 0;

    for(int i = 0; i < 4; i++){
        std::cout << "Round " << i+1 << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        Matrix C = A.mult(B);
        time += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
    }
    std::cout << "Time: " << time / 4 << "ms" << std::endl;
}