#include <iostream>
#include "Matrix.h"

int main(){
    std::vector<int> data = {1 , 2 , 3 , 4 , 5 ,
                             6 , 7 , 8 , 9 , 10,
                             11, 12, 13, 14, 15,
                             16, 17, 18, 19, 20,
                             21, 22, 23, 24, 25};
    Matrix<int> t(5,5,data);
    Matrix<int> test = t.cRange(1,3).rRange(1,3);
    return 0;
}