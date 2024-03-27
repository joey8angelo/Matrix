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
```

Matrix addition subtraction and multiplication.
```c++
D = A + B;
E = B - A;
F = A * B;
G = A.add(B);
H = A.sub(B);
I = B.mult(A);
```

Transpose.
```c++
A = A.T();
```