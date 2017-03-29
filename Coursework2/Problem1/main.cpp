/************************************************
 * in this script we'll solve
 * Boundary Value Problems of the form:
 *
 * u''(x) = f(x), on 0 < x < 1.
 * where u(0) = alpha, u(1) = beta
 *
 * using a Second Order Finite Difference Scheme.
 ***********************************************/
#include <iostream>
#include <cmath>
#include <string>
#include "Vector.hpp"
#include "Helpers.hpp"

static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/Coursework2/Problem1/OutputData/";
Helpers hlp;

double rhs_f(double x);

Vector true_solution(const Vector& mesh);

void solve(double (*rhs_f) (double x),
           double l_BC, double r_BC,
           const Vector& mesh);

void solve(double (*rhs_f) (double x),
           double l_BC, double r_BC,
           const Vector& mesh,
           Vector& result,
           bool save=false);

int main(){

    int n_points=9;
    //int n_points=200;
    double alpha = 0;   // Left Boundary Condition.
    double beta = 1;    // Right Boundary Condition.
    Vector mesh = hlp.linspace(0,1,n_points+1);

    // Right hand side function:
    double (*this_rhs) (double x);
    this_rhs = &rhs_f; // f(x) = 1

    // Algorithmic Implementation
    Vector sol(n_points+1);
    solve(this_rhs, alpha, beta, mesh, sol, true);

    Vector true_sol = true_solution(mesh);
    hlp.saveVectorToFile(true_sol, COMPLETE_PATH+"true.dat");


    // Code Verification
    std::ofstream this_file = hlp.openOutputFile(COMPLETE_PATH+"errors.dat");
    double mesh_size = 0.5;
    while(mesh_size > 1e-3){
        mesh_size *= 0.5;
        int this_n = 1/mesh_size;
        Vector this_mesh = hlp.linspace(0,1,this_n);
        Vector this_sol(this_n);
        solve(this_rhs, alpha, beta, this_mesh, this_sol);
        Vector this_true = true_solution(this_mesh);
        Vector difference = this_true-this_sol;
        double norm = difference.InfinityNorm();
        hlp.saveData(mesh_size, norm, this_file);
    }
    hlp.closeOutputFile(this_file);

    return 0;
}

// solve that saves the solution on a given vector.
void solve(double (*rhs_f) (double x),
           double l_BC, double r_BC,
           const Vector& mesh,
           Vector& result, bool save){

    int n = mesh.GetSize();

    if(n<=2)
        throw std::runtime_error("More than two mesh points should be given in order to compute a BVP.");

    if(result.GetSize() != n)
        throw std::runtime_error("Vector to save solution and mesh have differrent size.");
    int inner_n = n-2;

    // Assumes an equally spaced grid
    double h = mesh.Read(1) - mesh.Read(0);

    Vector diagonal(inner_n);
    Vector off_diag(inner_n-1);
    Vector rhs(inner_n);

    // Left Boundary condition on Right hand side
    rhs[0] = rhs_f( mesh.Read(1) ) - l_BC/pow(h,2);

    // Assemble inner tridiagonal "matrix"
    // and inner right hand side.
    for (int i = 0; i<inner_n-1; i++){
        diagonal[i] = -2.0/pow(h, 2);
        off_diag[i] = 1.0/pow(h,2);
        rhs[i+1] = rhs_f( mesh.Read(i+2) );
    }
    diagonal[inner_n-1] = -2.0/pow(h, 2);

    // Right Boundary condition on Right hand side
    rhs[inner_n-1] = rhs_f( mesh.Read(n-2) ) - r_BC/pow(h,2);

    // solve the system
    // rhs gets modified to contain the solution.
    hlp.solveTridiagonal(diagonal, off_diag, rhs);

    // Assemble the complete solution
    result[0] = l_BC;
    for (int i = 1; i<=inner_n; i++){
        result[i] = rhs[i-1];
    }
    result[n-1] = r_BC;

    // save everything on files.
    if(save){
        hlp.saveVectorToFile(mesh, COMPLETE_PATH+"mesh.dat");
        hlp.saveVectorToFile(result, COMPLETE_PATH+"computed.dat");
    }

}

double true_sol(double x){
    return x - sin(M_PI*x);
}

double rhs_f(double x){
    return pow(M_PI, 2)*sin(M_PI*x);
}

Vector true_solution(const Vector& xvec){
    int n = xvec.GetSize();
    Vector result(n);
    for (int i=0; i<n; i++){
        result[i] = true_sol(xvec.Read(i));
    }
    return result;
}
