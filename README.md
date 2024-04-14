# Matrix
A C++ implemenntation of a 2D Matrix that is focused on fast, IO efficient matrix multiplication and transpose. 
## Usage
The constructor for the Matrix class defines the number of columns and rows.
```c++
Matrix<int> A(2, 3);
```
Matrix A has 2 rows and 8 columns of integers.

A Matrix can also be defined with a default value.
```c++
Matrix<int> B(3,1,5);
```

Or a vector of data that aligns with the demensions.
```c++
std::vector<int> data = {1,2,3,4,5,6};
Matrix<int> C(2,3,data);
```

The shape of the matrix can be obtained with the shape() method, which returns a pair<int,int>.
A string representation of the shape is obtained with strShape() method.
```c++
std::cout << A.shape().first << ", " A.shape().second << std::endl; // prints 2,8
std::cout << A.strShape() << std::endl;                             // prints (2,8)
```

Access a position in a matrix using the at() method, which is zero-indexed.
```c++
A.at(0,0) = 1;
A.at(0,1) = 2;
A.at(0,2) = 3;
A.at(1,0) = 4;
A.at(1,1) = 5;
A.at(1,2) = 6;
/*
A = 1 2 3
    4 5 6
A == C
*/
```

The inner product of two matrices of compatible dimensions.
```c++
auto C = A.dot(B);
/*
A = 1 2 3   B = 5  C = 30
    4 5 6       5      75
                5
*/
```

Scalar addition subtraction and multiplication.
```c++
auto D = A + 5;
auto E = 5 - B;
auto F = A * 5;
auto G = A.add(5);
auto H = B.sub(5);
auto I = A.mult(5);
/*
D = 6 7  8   E = 0  F = 5  10 15
    9 10 11      0      20 25 30
                 0
*/
```

Matrix addition subtraction and multiplication.
```c++
D = A + A;
E = B - B;
F = A * A;
G = A.add(A);
H = B.sub(B);
I = A.mult(A);
/*
D = 2 4  6   E = 0   F = 1  4  9
    8 10 12      0       16 25 36
                 0
*/
```

Transpose.
```c++
A = A.T();
/*
A = 1 4
    2 5
    3 6
*/
```

And other useful functions, reshape, range access, log, exp, max, argMax.

## Compiling
Compile with -O3 -march=native on widnows, or use the makeFile, which compiles a main.cpp file by default. 
