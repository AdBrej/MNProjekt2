#include <iostream>
#include "Eigen/Dense"
#include "constValue.h"
#include <math.h>
#include <iomanip>
#include <chrono>

struct Matrix{
    double**Matrix;
    int sizeRow;
    int sizeCol;
};

struct Result{
    Matrix M;
    int iteration;
    double norm;
    long long duration;

};
using Eigen::MatrixXd;


Matrix createMatrix(int sizeRow, int sizeCol, double inValue = 0){
    double** A = 0;
    A = new double*[sizeRow];
    for (int h = 0; h < sizeRow; h++)
    {
        A[h] = new double[sizeCol];
        for (int w = 0; w < sizeCol; ++w) {
            A[h][w] = inValue;
        }
    }
    Matrix M;
    M.sizeRow = sizeRow;
    M.sizeCol = sizeCol;
    M.Matrix = A;
    return M;
}




void displayMatrix(Matrix &M){
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    for (int row = 0; row < sizeRow; ++row) {
        std::cout << std::setw(2) << std::right << row << " | ";
        for (int col = 0; col < sizeCol; ++col) {
            std::cout << std::setw(2) << std::right << M.Matrix[row][col] << " ";
        }
        std::cout << std::endl;
    }
}

Matrix addMatrixNxN(Matrix &M, Matrix &M2, int modyfikator = 1){
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    Matrix out = createMatrix(sizeRow, sizeCol);
    for (int row = 0; row < sizeRow; ++row) {
        for (int col = 0; col < sizeCol; ++col) {
            out.Matrix[row][col] = M.Matrix[row][col] + (M2.Matrix[row][col] * modyfikator);
        }
    }
    return out;
}






Matrix DMatrix(Matrix &M){
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    double** D = 0;
    D = new double*[sizeRow];
    for (int h = 0; h < sizeRow; h++)
    {
        D[h] = new double[sizeCol];
        for (int w = 0; w < sizeCol; ++w) {
            if(h == w){
                D[h][w] = M.Matrix[h][w];
            }else{

                D[h][w] = 0;
            }
        }
    }
    Matrix out;
    out.Matrix = D;
    out.sizeCol = sizeCol;
    out.sizeRow = sizeRow;
    return out;
}



Matrix LMatrix(Matrix &M){
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    double** L = 0;
    L = new double*[sizeRow];
    for (int h = 0; h < sizeRow; h++)
    {
        L[h] = new double[sizeCol];
        for (int w = 0; w < sizeCol; ++w) {
            if(w < h and M.Matrix[h][w] != 0){
                L[h][w] = M.Matrix[h][w] * (-1);
            }else{

                L[h][w] = 0;
            }
        }
    }
    Matrix out;
    out.Matrix = L;
    out.sizeCol = sizeCol;
    out.sizeRow = sizeRow;
    return out;
}

Matrix UMatrix(Matrix &M){
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    double ** U = 0;
    U = new double*[sizeRow];
    for (int h = 0; h < sizeRow; h++)
    {
        U[h] = new double[sizeCol];
        for (int w = 0; w < sizeCol; ++w) {
            if(w > h and M.Matrix[h][w] != 0){
                U[h][w] = M.Matrix[h][w] * (-1);
            }else{

                U[h][w] = 0;
            }
        }
    }
    Matrix out;
    out.Matrix = U;
    out.sizeRow = sizeRow;
    out.sizeCol = sizeCol;
    return out;
}


Matrix createMyMatrixNxN(int size, double a1, double a2, double a3){

    Matrix out = createMatrix(size, size);

    for (int i = 0; i < size ; ++i) {
        out.Matrix[i][i] = a1;
        if (i + 1 < size) {
            out.Matrix[i][i+1] = a2;
            out.Matrix[i+1][ i] = a2;
        }
        if (i + 2 < size) {
            out.Matrix[i][i+2] = a3;
            out.Matrix[i+2][ i] = a3;
        }
    }
    return out;

}
void createDiagonalNxN(Matrix &M, double value = 1){
    int sizeRow = M.sizeRow;

    for (int row = 0; row < sizeRow; ++row) {
        M.Matrix[row][row] = value;
    }
}
Matrix copyMatrix(Matrix &M){
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    Matrix out = createMatrix(sizeRow, sizeCol);
    for (int row = 0; row < sizeRow; ++row) {
        for (int col = 0; col < sizeCol; ++col) {
            out.Matrix[row][col] = M.Matrix[row][col];
        }
    }
    return out;
}
Matrix inverseMatrix(Matrix &A){
    int sizeRow = A.sizeRow;
    int sizeCol = A.sizeCol;
    Matrix I = createMatrix(sizeRow, sizeCol);
    createDiagonalNxN(I);
    Matrix tmp = copyMatrix(A);
    double a, b; // współczynniki przy macierzach sprawdzone przed zmianamina wierszu
    for (int row = 0; row < sizeRow; ++row) {
        a = tmp.Matrix[row][row];
        for (int col = 0; col < sizeCol; ++col) {
            if(a== 0){
                break;
            }else{
                double tmp1 = tmp.Matrix[row][col] / a;
                tmp.Matrix[row][col] = tmp1;
                double tmp2 = I.Matrix[row][col] / a;
                I.Matrix[row][col] = tmp2;
            }
        }
        for (int row2 = 0; row2 < sizeRow; ++row2) {
            if(row2 != row){
                b = tmp.Matrix[row2][row];
                for (int col = 0; col < sizeCol; ++col) {
                    double tmp1 = tmp.Matrix[row2][col] - (tmp.Matrix[row][col] * b);
                    tmp.Matrix[row2][col] = tmp1;
                    double tmp2 = I.Matrix[row2][col] - (I.Matrix[row][col] * b);
                    I.Matrix[row2][col] = tmp2;
                }
            }
        }
    }

    return I;
}

Matrix multiplyMatrix(Matrix &M, Matrix &X){
    int MsizeRow = M.sizeRow;
    int MsizeCol = M.sizeCol;
    int XsizeRow = X.sizeRow;
    int XsizeCol = X.sizeCol;
    Matrix out = createMatrix(XsizeRow, XsizeCol);
    for (int rowM = 0; rowM < MsizeRow; ++rowM) {
        for (int colX = 0; colX < XsizeCol; ++colX) {
            double sum = 0;
            for (int i = 0; i < MsizeCol; ++i) {
                double tmp = M.Matrix[rowM][i] * X.Matrix[i][colX];
                sum += tmp;
            }
            out.Matrix[rowM][colX] = sum;
        }
    }
    return out;
}



Matrix createMyVectorN(int size, double f){
    Matrix out = createMatrix(size, 1);
    for (int i = 0; i < size ; ++i) {
        out.Matrix[i][0] = std::sin((i*(f + 1)/50));
    }
    return out;
}


Matrix createVectorN(int size, double value){
    Matrix out = createMatrix(size, 1);
    for (int i = 0; i < size ; ++i) {
        out.Matrix[i][0] = value;
    }
    return out;
}


double norm(Matrix &M){
    double sum = 0;
    int sizeRow = M.sizeRow;
    int sizeCol = M.sizeCol;
    for (int row = 0; row < sizeRow; ++row) {
        for (int col = 0; col < sizeCol; ++col) {
            sum += M.Matrix[row][col] * M.Matrix[row][col];
        }
    }
    return sqrt(sum);
}

Result GaussS(Matrix &A, Matrix &b, double eps){
    long long ms = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    Matrix Xin = createMatrix(b.sizeRow, b.sizeCol);
    int iterationCount = 0;

    Matrix D = DMatrix(A);
    Matrix L = LMatrix(A);
    Matrix U = UMatrix(A);
//    Matrix DL = addMatrixNxN(D, L, -1);
    Matrix DL = addMatrixNxN(D, L, mConst::ODEJMOWANIE);

    Matrix iDL = inverseMatrix(DL);
    Matrix T = multiplyMatrix(iDL, U);
    Matrix F = multiplyMatrix(iDL, b);
//    displayMatrix(A);
//    displayMatrix(b);
//    displayMatrix(D);
//    displayMatrix(L);
//    displayMatrix(U);
//    displayMatrix(DL);
//    displayMatrix(iDL);
//    displayMatrix(T);
//    displayMatrix(F);

    Matrix Xout;
    Matrix TX;
    Matrix XoutXin;
    double nor;

    do {
        TX = multiplyMatrix(T, Xin);
        Xout = addMatrixNxN(TX, F);
//        XoutXin = addMatrixNxN(Xout, Xin, -1);
        XoutXin = addMatrixNxN(Xout, Xin, mConst::ODEJMOWANIE);
        iterationCount += 1;
        Xin = Xout;
        nor = norm(XoutXin);
    } while (eps < nor);
    long long ms2 = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    Result result;
    result.duration = ms2 - ms;
    result.M = Xout;
    result.iteration = iterationCount;
    result.norm = nor;
    return result;
}

Result Jacobi(Matrix &A, Matrix &b, double eps) {
    long long ms = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    Matrix Xin = createMatrix(b.sizeRow, b.sizeCol);
    int iterationCount = 0;

    Matrix D = DMatrix(A);
    Matrix L = LMatrix(A);
    Matrix U = UMatrix(A);
    Matrix LU = addMatrixNxN(L, U);
    Matrix iD = inverseMatrix(D);
    Matrix T = multiplyMatrix(iD, LU);
    Matrix F = multiplyMatrix(iD, b);

    Matrix Xout;
    Matrix TX;
    Matrix XoutXin;
    double nor;

    do {
        TX = multiplyMatrix(T, Xin);
        Xout = addMatrixNxN(TX, F);
//        XoutXin = addMatrixNxN(Xout, Xin, -1);
        XoutXin = addMatrixNxN(Xout, Xin, mConst::ODEJMOWANIE);
        iterationCount += 1;
        Xin = Xout;
        nor = norm(XoutXin);
    } while (eps < nor);
    long long ms2 = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    Result result;
    result.duration = ms2 - ms;
    std::cout << ms2 <<'\n'<< ms<<'\n';
    result.M = Xout;
    result.iteration = iterationCount;
    result.norm = nor;
    return result;
}


Result Gaussa(Matrix &A, Matrix &B){
    long long ms = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    Matrix I = copyMatrix(B);
//    createDiagonalNxN(I);
    int iteration = 0;
    Matrix tmp = copyMatrix(A);
    double a, b; // współczynniki przy macierzach sprawdzone przed zmianamina wierszu
    for (int row = 0; row < tmp.sizeRow; ++row) {
        a = tmp.Matrix[row][row];
        iteration +=1;
        std::cout << iteration <<"\n";
        for (int col = 0; col < tmp.sizeCol; ++col) {
            if(a== 0){
                break;
            }else{
                double tmp1 = tmp.Matrix[row][col] / a;
                tmp.Matrix[row][col] = tmp1;
//                double tmp2 = I.Matrix[row][col] / a;
//                I.Matrix[row][col] = tmp2;
            }
        }
        if(a != 0){
            double tmp2 = I.Matrix[row][0] / a;
            I.Matrix[row][0] = tmp2;
        }
        for (int row2 = 0; row2 < tmp.sizeRow; ++row2) {
            if(row2 != row){
                b = tmp.Matrix[row2][row];
                for (int col = 0; col < tmp.sizeCol; ++col) {
                    double tmp1 = tmp.Matrix[row2][col] - (tmp.Matrix[row][col] * b);
                    tmp.Matrix[row2][col] = tmp1;
//                    double tmp2 = I.Matrix[row2][0] - (I.Matrix[row][0] * b);
//                    I.Matrix[row2][0] = tmp2;
                }
                double tmp2 = I.Matrix[row2][0] - (I.Matrix[row][0] * b);
                I.Matrix[row2][0] = tmp2;
            }
        }
    }

    long long ms2 = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    Result result;
    result.duration = ms2 - ms;
    result.M = I;
    result.iteration = 0;
//    result.norm = norm(addMatrixNxN(multiplyMatrix(A, I), B, -1));
    Matrix multiAandI = multiplyMatrix(A, I);
    Matrix addMforNorm = addMatrixNxN(multiAandI, B, mConst::ODEJMOWANIE);
    result.norm = norm(addMforNorm);
//    displayMatrix(tmp);
//    displayMatrix(I);
    return result;

}


int main() {
//    std::cout.precision(4);
    std::cout.width( 10 );

//zad1 A
//    MatrixXd A;
//    A = createMatrixNxN(mConst::N,  5 + mConst::e,-1,-1);
//    Eigen::VectorXd b;
//    b = createVectorN(30, mConst::f);



//    testowe mniejsze rozmiary by było widać co sie dzieje przy debugowaniu

    Matrix A = createMyMatrixNxN(mConst::tmpsize, (double) 5 + mConst::e, (double) -1, (double) -1);

    Matrix b = createMyVectorN(mConst::tmpsize, (double) mConst::f);

//    displayMatrix(A);
//    displayMatrix(b);
/*

    MatrixXd z(3,3);
    double g[3][3] = {-1,-1,-1,2,1,4,1,1,2};

    double ** y = 0;
    y = new double*[3];
    for (int h = 0; h < 3; h++)
    {
        y[h] = new double[3];
        for (int w = 0; w < 3; ++w) {
            y[h][w] = g[h][w];
        }
    }
    Matrix F;
    F.Matrix = y;
    F.sizeCol = 3;
    F.sizeRow = 3;
    z << -1, -1, -1,
            2, 1, 4,
            1, 1, 2;

    double g2[2][3] = {1,0,2,-1,3,1};
    double g3[3][2] = {3, 1, 2, 1, 1, 0};

    double ** y2 = 0;
    y2 = new double*[2];
    for (int h = 0; h < 2; h++)
    {
        y2[h] = new double[3];
        for (int w = 0; w < 3; ++w) {
            y2[h][w] = g2[h][w];
        }
    }
    Matrix F2;
    F2.Matrix = y2;
    F2.sizeCol = 3;
    F2.sizeRow = 2;

    double ** y3 = 0;
    y3 = new double*[3];
    for (int h = 0; h < 3; h++)
    {
        y3[h] = new double[2];
        for (int w = 0; w < 2; ++w) {
            y3[h][w] = g3[h][w];
        }
    }
    Matrix F3;
    F3.Matrix = y3;
    F3.sizeCol = 2;
    F3.sizeRow = 3;
    std::cout << z << std::endl;
//    std::cout.precision(10);
    std::cout << z.inverse() << std::endl;
    //displayMatrixT(3, DMatrixT(3, y));
    //displayMatrixT(F.sizeRow, UMatrix(F.sizeRow, F.Matrix));

    Matrix DDD = inverseMatrix(F);
    //displayMatrixT(DDD.sizeCol, DDD.Matrix);
//    displayMatrixT(3, y);
//    inversMatrixNxN(3, y);
//    displayMatrixT(3, inversMatrixNxN(3, y));

    std::cout << "--------------------------" << std::endl;
    Matrix F4 = multiplyMatrix(F2, F3);
    displayMatrixT(2, multiplyMatrix(F2, F3).Matrix);
//    Matrix s = addMatrixNxN(F, F, -1);
    Matrix s = addMatrixNxN(F, F, mConst::ODEJMOWANIE);
    displayMatrixT(s.sizeRow, s.Matrix);
    Matrix AA = createMatrix(4, 4);
    Matrix BB = createMatrix(4, 1);
    double g4[4][4] = {10, -1, 2, 0,
                       -1, 11, -1, 3,
                        2, -1, 10, -1,
                        0, 3, -1, 8};

    double ** y4 = 0;
    y4 = new double*[4];
    for (int h = 0; h < 4; h++)
    {
        y4[h] = new double[4];
        for (int w = 0; w < 4; ++w) {
            y4[h][w] = g4[h][w];
        }
    }
    AA.Matrix = y4;
    BB.Matrix[0][0] = 6;
    BB.Matrix[1][0] = 25;
    BB.Matrix[2][0] = -11;
    BB.Matrix[3][0] = 15;

    double eps= 1e-9;
    Result RRR = GaussS(AA, BB, eps);
    Result RRR2 = Jacobi(AA, BB, eps);
    Result RRR3 = Gaussa(AA, BB);
    displayMatrix(RRR.M);
    std::cout << RRR.iteration;
    displayMatrix(RRR2.M);
    std::cout << RRR2.iteration;
    displayMatrix(RRR3.M);
    std::cout << RRR3.iteration;
    std::cout << "--------------------------" << std::endl;
    Gaussa(AA, BB);
    inverseMatrix(AA);
    long long ms = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << ms;
*/

    Result R = Jacobi(A, b, mConst::eps);
    std::cout<<"----------------------------------------------------" <<std::endl;
    std::cout<<R.iteration<<"\n"<<R.duration<<"\n";
    displayMatrix(R.M);
    Result R2 = GaussS(A, b, mConst::eps);
    std::cout<<"----------------------------------------------------" <<std::endl;
    std::cout << R2.iteration << "\n" << R2.duration << "\n";
    displayMatrix(R2.M);
    Result R3 = Gaussa(A, b);
    std::cout<<"----------------------------------------------------" <<std::endl;
    std::cout << R3.iteration << "\n" << R3.duration << "\n";
    displayMatrix(R3.M);
//    ----------KONIEC TESTÓW--------------------Matrix AA = createMatrix(4, 4);
    Matrix BB = createMatrix(4, 1);
    BB.Matrix[0][0] = 13.15;
    BB.Matrix[1][0] = 49.84;
    BB.Matrix[2][0] = -14.08;
    BB.Matrix[3][0] = -46.51;
    Matrix AA = createMatrix(4,4);
    double g4[4][4] = {1.2, 2.6, -0.1, 1.5,
                       4.5, 9.8, -0.4, 5.7,
                       0.1, -0.1, -0.3, -3.5,
                       4.5, -5.2, 4.2, -3.4};

    double ** y4 = 0;
    y4 = new double*[4];
    for (int h = 0; h < 4; h++)
    {
        y4[h] = new double[4];
        for (int w = 0; w < 4; ++w) {
            y4[h][w] = g4[h][w];
        }
    }
    AA.Matrix = y4;
    Result R4 = Gaussa(AA, BB);
    std::cout<<"----------------------------------------------------" <<std::endl;
    std::cout << R4.iteration << "\n" << R4.duration << "\n";
    displayMatrix(R4.M);
    return 0;
}