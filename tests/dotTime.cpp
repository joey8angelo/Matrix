#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include "../Matrix.h"

int main(){
    std::vector<int> Xtime;
    std::vector<int> Ytime;

    for(int i = 5; i <= 2000; i+=10){
        Matrix<double> A(i,i);
        Matrix<double> B(i,i);
        srand(0);
        for(int a = 0; a < i; a++){
            for(int b = 0; b < i; b++){
                A.at(a,b) = rand() % 10;
                B.at(a,b) = rand() % 10;
            }
        }
        Xtime.push_back(i);
        int time = 0;
        std::cout << "Testing size: " << i << std::endl;
        for(int r = 0; r < 3; r++){
            auto start = std::chrono::high_resolution_clock::now();
            Matrix C = A.mult(B);
            time += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        }
        std::cout << "Time: " << time / 3 << std::endl;
        Ytime.push_back(time / 3);
    }

    std::cout << "Xtime = [";
    for(std::size_t i = 0; i < Xtime.size(); i++){
        if(i)
            std::cout << ", " << Xtime[i] << " ";
        else
            std::cout << Xtime[i] << " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Ytime = [";
    for(std::size_t i = 0; i < Ytime.size(); i++){
        if(i)
            std::cout << ", " << Ytime[i] << " ";
        else
            std::cout << Ytime[i] << " ";
    }
    std::cout << "]" << std::endl;
}