#include <iostream>
#include <vector>
#include <random>
#include "../Matrix.h"

std::vector<std::vector<double>> dot(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t k = 0; k < A[0].size(); k++){
            for(std::size_t j = 0; j < B[0].size(); j++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>>& A){
    std::vector<std::vector<double>> C(A[0].size(), std::vector<double>(A.size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[j][i] = A[i][j];
        }
    }
    return C;
}

std::vector<std::vector<double>> add(std::vector<std::vector<double>>& A, double B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] + B;
        }
    }
    return C;
}

std::vector<std::vector<double>> add(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<double>> sub(std::vector<std::vector<double>>& A, double B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] - B;
        }
    }
    return C;
}

std::vector<std::vector<double>> sub(double A, std::vector<std::vector<double>>& B){
    std::vector<std::vector<double>> C(B.size(), std::vector<double>(B[0].size(), 0));
    for(std::size_t i = 0; i < B.size(); i++){
        for(std::size_t j = 0; j < B[0].size(); j++){
            C[i][j] = A - B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<double>> sub(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<double>> mult(std::vector<std::vector<double>>& A, double B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] * B;
        }
    }
    return C;
}

std::vector<std::vector<double>> mult(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B){
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] * B[i][j];
        }
    }
    return C;
}

bool comp(Matrix<double>& A, std::vector<std::vector<double>>& Av, Matrix<double>& B, std::vector<std::vector<double>>& Bv, int op, int ops){
    for(int a = 0; a < A.shape().first; a++){
        for(int b = 0; b < A.shape().second; b++){
            if(A.at(a,b) != Av[a][b]){
                std::cerr << "Error in operation " << op << " after " << ops << " operations at A[" << a << "][" << b << "]: " << A.at(a,b) << " != " << Bv[a][b] << std::endl;
                return true;
            }
            if(B.at(a,b) != Bv[a][b]){
                std::cerr << "Error in operation " << op << " after " << ops << " operations at B[" << a << "][" << b << "]: " << B.at(a,b) << " != " << Bv[a][b] << std::endl;
                return true;
            }
        }
    }
    return false;
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
    srand(t);
    std::cout << "Testing operations " << tests << " times with seed " << t << std::endl;
    for(long int i = 0; i < tests; i++){
        int N = rand() % 3000 + 1;
        int M = rand() % 3000 + 1;

        Matrix<double> A(N, M);
        Matrix<double> B(N, M);

        std::vector<std::vector<double>> Av(N, std::vector<double>(M));
        std::vector<std::vector<double>> Bv(N, std::vector<double>(M));

        for(int a = 0; a < N; a++){
            for(int b = 0; b < M; b++){
                Av[a][b] = rand() % 10;
                A.at(a,b) = Av[a][b];
            }
        }
        for(int a = 0; a < N; a++){
            for(int b = 0; b < M; b++){
                Bv[a][b] = rand() % 10;
                B.at(a,b) = Bv[a][b];
            }
        }

        for(int ops = 0; ops < 15; ops++){
            int op = rand() % 25;
            switch(op){
                case 0:
                    A = A + B;
                    Av = add(Av, Bv);
                    break;
                case 1:
                    A = B + A;
                    Av = add(Bv, Av);
                    break;
                case 2:
                    B = A + B;
                    Bv = add(Av, Bv);
                    break;
                case 3:
                    B = B + A;
                    Bv = add(Bv, Av);
                    break;
                case 4:
                    A = A - B;
                    Av = sub(Av, Bv);
                    break;
                case 5:
                    A = B - A;
                    Av = sub(Bv, Av);
                    break;
                case 6:
                    B = A - B;
                    Bv = sub(Av, Bv);
                    break;
                case 7:
                    B = B - A;
                    Bv = sub(Bv, Av);
                    break;
                case 8:
                    A = A * B;
                    Av = mult(Av, Bv);
                    break;
                case 9:
                    A = B * A;
                    Av = mult(Bv, Av);
                    break;
                case 10:
                    B = A * B;
                    Bv = mult(Av, Bv);
                    break;
                case 11:
                    B = B * A;
                    Bv = mult(Bv, Av);
                    break;
                case 12:
                    A = A + 1;
                    Av = add(Av, 1);
                    break;
                case 13:
                    A = 1 + A;
                    Av = add(Av, 1);
                    break;
                case 14:
                    B = B + 1;
                    Bv = add(Bv, 1);
                    break;
                case 15:
                    B = 1 + B;
                    Bv = add(Bv, 1);
                    break;
                case 16:
                    A = A - 1;
                    Av = sub(Av, 1);
                    break;
                case 17:
                    A = 1 - A;
                    Av = sub(1, Av);
                    break;
                case 18:
                    B = B - 1;
                    Bv = sub(Bv, 1);
                    break;
                case 19:
                    B = 1 - B;
                    Bv = sub(1, Bv);
                    break;
                case 20:
                    A = A * 2;
                    Av = mult(Av, 2);
                    break;
                case 21:
                    A = 2 * A;
                    Av = mult(Av, 2);
                    break;
                case 22:
                    B = B * 2;
                    Bv = mult(Bv, 2);
                    break;
                case 23:
                    B = 2 * B;
                    Bv = mult(Bv, 2);
                    break;
                case 24:
                    A = A.T();
                    Av = transpose(Av);
                    B = B.T();
                    Bv = transpose(Bv);
                    break;
            }
            if(comp(A, Av, B, Bv, op, ops)){
                return 1;
            }
        }
        std::cout << "Test " << i << " passed" << std::endl;
    }
}