#include <iostream>
#include <vector>
#include <random>
#include "Matrix.h"

std::vector<std::vector<float>> add(std::vector<std::vector<float>>& A, float B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] + B;
        }
    }
    return C;
}

std::vector<std::vector<float>> add(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<float>> sub(std::vector<std::vector<float>>& A, float B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] - B;
        }
    }
    return C;
}

std::vector<std::vector<float>> sub(float A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(B.size(), std::vector<float>(B[0].size(), 0));
    for(std::size_t i = 0; i < B.size(); i++){
        for(std::size_t j = 0; j < B[0].size(); j++){
            C[i][j] = A - B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<float>> sub(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<float>> mult(std::vector<std::vector<float>>& A, float B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] * B;
        }
    }
    return C;
}

std::vector<std::vector<float>> mult(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] * B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<float>> div(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] / B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<float>> div(std::vector<std::vector<float>>& A, float B){
    std::vector<std::vector<float>> C(A.size(), std::vector<float>(A[0].size(), 0));
    for(std::size_t i = 0; i < A.size(); i++){
        for(std::size_t j = 0; j < A[0].size(); j++){
            C[i][j] = A[i][j] / B;
        }
    }
    return C;
}

std::vector<std::vector<float>> div(float A, std::vector<std::vector<float>>& B){
    std::vector<std::vector<float>> C(B.size(), std::vector<float>(B[0].size(), 0));
    for(std::size_t i = 0; i < B.size(); i++){
        for(std::size_t j = 0; j < B[0].size(); j++){
            C[i][j] = A / B[i][j];
        }
    }
    return C;
}

bool comp(Matrix<float>& A, std::vector<std::vector<float>>& Av, Matrix<float>& B, std::vector<std::vector<float>>& Bv, int op, int ops){
    for(int a = 0; a < A.shape().first; a++){
        for(int b = 0; b < A.shape().second; b++){
            // NaN == NaN is false, so check for NaNs separately, just to ensure that the Matrix is behaving correctly
            bool aNaN = std::isnan(A.at(a,b));
            bool avNaN = std::isnan(Av[a][b]);
            if(aNaN || avNaN){
                if(aNaN != avNaN){
                    std::cerr << "Error in operation " << op << " after " << ops << " operations at A[" << a << "][" << b << "]: " << A.at(a,b) << " != " << Av[a][b] << std::endl;
                    return true;
                }
            }
            else if(A.at(a,b) != Av[a][b]){
                std::cerr << "Error in operation " << op << " after " << ops << " operations at A[" << a << "][" << b << "]: " << A.at(a,b) << " != " << Av[a][b] << std::endl;
                return true;
            }
            bool bNaN = std::isnan(B.at(a,b));
            bool bvNaN = std::isnan(Bv[a][b]);
            if(bNaN || bvNaN){
                if(bNaN != bvNaN){
                    std::cerr << "Error in operation " << op << " after " << ops << " operations at B[" << a << "][" << b << "]: " << B.at(a,b) << " != " << Bv[a][b] << std::endl;
                    return true;
                }
            }
            else if(B.at(a,b) != Bv[a][b]){
                std::cerr << "Error in operation " << op << " after " << ops << " operations at A[" << a << "][" << b << "]: " << B.at(a,b) << " != " << Bv[a][b] << std::endl;
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
        int M = rand() % 5 + 1;
        int N = rand() % 5 + 1;

        Matrix<float> A(M, N);
        Matrix<float> B(M, N);

        std::vector<std::vector<float>> Av(M, std::vector<float>(N));
        std::vector<std::vector<float>> Bv(M, std::vector<float>(N));

        for(int a = 0; a < M; a++){
            for(int b = 0; b < N; b++){
                Av[a][b] = rand() % 1000 - 500;
                A.at(a,b) = Av[a][b];
            }
        }
        for(int a = 0; a < M; a++){
            for(int b = 0; b < N; b++){
                Bv[a][b] = rand() % 1000 - 500;
                B.at(a,b) = Bv[a][b];
            }
        }

        for(int ops = 0; ops < 15; ops++){
            int op = rand() % 32;
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
                    A = A / B;
                    Av = div(Av, Bv);
                    break;
                case 13:
                    A = B / A;
                    Av = div(Bv, Av);
                    break;
                case 14:
                    B = A / B;
                    Bv = div(Av, Bv);
                    break;
                case 15:
                    B = B / A;
                    Bv = div(Bv, Av);
                    break;
                case 16:
                    A = A + 1;
                    Av = add(Av, 1);
                    break;
                case 17:
                    A = 1 + A;
                    Av = add(Av, 1);
                    break;
                case 18:
                    B = B + 1;
                    Bv = add(Bv, 1);
                    break;
                case 19:
                    B = 1 + B;
                    Bv = add(Bv, 1);
                    break;
                case 20:
                    A = A - 1;
                    Av = sub(Av, 1);
                    break;
                case 21:
                    A = 1 - A;
                    Av = sub(1, Av);
                    break;
                case 22:
                    B = B - 1;
                    Bv = sub(Bv, 1);
                    break;
                case 23:
                    B = 1 - B;
                    Bv = sub(1, Bv);
                    break;
                case 24:
                    A = A * 2;
                    Av = mult(Av, 2);
                    break;
                case 25:
                    A = 2 * A;
                    Av = mult(Av, 2);
                    break;
                case 26:
                    B = B * 2;
                    Bv = mult(Bv, 2);
                    break;
                case 27:
                    B = 2 * B;
                    Bv = mult(Bv, 2);
                    break;
                case 28:
                    A = A / 2;
                    Av = div(Av, 2);
                    break;
                case 29:
                    A = 2 / A;
                    Av = div(2, Av);
                    break;
                case 30:
                    B = B / 2;
                    Bv = div(Bv, 2);
                    break;
                case 31:
                    B = 2 / B;
                    Bv = div(2, Bv);
                    break;
            }
            if(comp(A, Av, B, Bv, op, ops)){
                return 1;
            }
        }
        std::cout << "Test " << i << " passed" << std::endl;
    }
}