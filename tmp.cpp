/*
//
// Created by Brej on 09.04.2017.
//

#include <iostream>
#include "Eigen/Dense"
#include "constValue.h"
#include <math.h>
using Eigen::MatrixXd;

template <typename T>
MatrixXd createMatrixNxN(int size,T a1,T a2,T a3){
    MatrixXd A;
    A = MatrixXd::Zero(size,size);
    for (int i = 0; i < size ; ++i) {
        A(i,i) = a1;
        if (i + 1 < size) {
            A(i,i+1) = a2;
            A(i+1, i) = a2;
        }
        if (i + 2 < size) {
            A(i,i+2) = a3;
            A(i+2, i) = a3;
        }
    }
    return A;
}

template <typename T>
Eigen::VectorXd createVectorN(int size, T f){
    Eigen::VectorXd b(size);
    for (int i = 0; i < size ; ++i) {
        b(i) = std::sin(i*(f + 1)/50);
    }
    return b;
}


int main() {
    std::cout << "Hello, World!" << std::endl;

//zad1 A
//    MatrixXd A;
//    A = createMatrixNxN(mConst::N,  5 + mConst::e,-1,-1);
//    Eigen::VectorXd b;
//    b = createVectorN(30, mConst::f);



//    testowe mniejsze rozmiary by było widać co sie dzieje przy debugowaniu
    MatrixXd A;
    A = createMatrixNxN(40,  5 + mConst::e,-1,-1);
    Eigen::VectorXd b;
    b = createVectorN(40, mConst::f);


    std::cout << A << std::endl;
    std::cout << b << std::endl;

    MatrixXd z(3,3);

    z << 1, 2, 1,
            2, 1, 0,
            -1, 1, 2;

    std::cout << z << std::endl;

    std::cout << z.inverse() << std::endl;
//    ----------KONIEC TESTÓW--------------------
    return 0;
}*/
