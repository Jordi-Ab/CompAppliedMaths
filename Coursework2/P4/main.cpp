#include <iostream>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "ProjectedSOR.hpp"

Vector createB(int size);

int main(){
    int n = 15; // Size of the system.
    int l = 1; // First entry of lower diagonal.
    int u = 1; // First entry of upper diagonal.
    int h = 2/(n + 1); // Mesh Size.
    int m = -2*(1+2*(h*h));

    double omega = 1.527893;

    Matrix A2(n,n);
    A2.tridiagonal(l,m,u);

    //std::cout << "Matrix:" << std::endl;
    //std::cout << A2 << std::endl;

    Vector b = createB(n);
    //std::cout << "b:" << std::endl;
    //std::cout << b << std::endl;

    Vector init_x(n);

    Matrix A(4,4);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 4; A(1,4) = 0;
    A(2,1) = 0; A(2,2) = 1; A(2,3) = 0; A(2,4) = 2;
    A(3,1) = 7; A(3,2) = 0; A(3,3) = 1; A(3,4) = 0;
    A(4,1) = 4; A(4,2) = 2; A(4,3) = 5; A(4,4) = 0;

    Vector x(n);

    ProjectedSOR lin_solver(A2, b, omega, init_x);
    lin_solver.SorSolve(x);
    //lin_solver.GsSolve(x);

    std::cout << "x:" << std::endl;
    std::cout << x << std::endl;

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

