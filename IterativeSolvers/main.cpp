#include <iostream>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "IterativeSolver.hpp"

Vector createB(int size);
void something();

int main(){
    int n = 5; // Size of the system.
    int l = 1; // First entry of lower diagonal.
    int u = 1; // First entry of upper diagonal.
    int m = -2;

    double omega = 1.527893;

    Matrix A2(n,n);
    A2.tridiagonal(l,m,u);

    std::cout << A2 << std::endl;

    Vector ld(n-1);
    ld.fillWith(l);
    std::cout << ld << std::endl;

    Vector d(n);
    d.fillWith(m);
    std::cout << d << std::endl;

    Vector ud(n-1);
    ud.fillWith(u);
    std::cout << ud << std::endl;

    //std::cout << "Matrix:" << std::endl;
    //std::cout << A2 << std::endl;

    Vector b = createB(n);
    //std::cout << "b:" << std::endl;
    //std::cout << b << std::endl;
    std::cout << b << std::endl;

    Vector init_x(n); // Zeros

    Vector x_sor1(n); // Will hold the result
    Vector x_sor2(n); // Will hold the result

    IterativeSolver lin_solverM(A2, b, init_x);
    IterativeSolver lin_solverT(ld, d, ud, b, init_x);

    lin_solverM.solveSOR(omega, x_sor1);
    lin_solverT.solveSOR(omega, x_sor2);

    std::cout << "Answer1:" << std::endl;
    std::cout << x_sor1 << std::endl;

    std::cout << "Answer2:" << std::endl;
    std::cout << x_sor2 << std::endl;

    return 0;
}

Vector createB(int size){
    Vector b(size);
    for(int entry=0; entry<size; entry++){
        if (entry%2 == 0){
            b[entry] = -1;
        }else{
            b[entry] = 0;
        }
    }
    return b;
}


