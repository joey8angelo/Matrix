#include <iostream>
#include <vector>
#include "Matrix.h"

int main(){
    std::vector<int> d = {1, 5, 3, 4, 2, 6};
    Matrix<int> m(2,3,d);
    auto t = m.max(0);
    return 0;
}