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

// Complete path of the folder where you want the solution to be saved.
//static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/CompAppliedMaths_git/Coursework2/P1/OutputData/";
static const std::string COMPLETE_PATH = "OutputData/";
Helpers hlp; // A set of useful functions.

/*
 * Function: rhs_f1
 * ---------------------------------------
 * Right hand side of the second order BVP.
 * U''(x) = f(x).
 *
 * Function on Case 1 of Coursework.
*/
double rhs_f1(double x);

/*
 * Function: rhs_f2
 * ---------------------------------------
 * Right hand side of the second order BVP.
 * U''(x) = f(x).
 *
 * Function on Case 2 of Coursework.
*/
double rhs_f2(double x);

/*
 * Function: true_solution1
 * ---------------------------------------------
 * True solution of the BVP for Case 1 on Coursework.
 * Note: Is a vector valued function.
 * Receives:
 *  mesh    -> Vector of points where to evaluate.
 * Returns:
 *  Vector of solution.
*/
Vector true_solution1(const Vector& mesh);

/*
 * Function: true_solution2
 * ---------------------------------------------
 * True solution of the BVP for Case 2 on Coursework.
 * Note: Is a vector valued function.
 * Receives:
 *  mesh    -> Vector of points where to evaluate.
 * Returns:
 *  Vector of solution.
*/
Vector true_solution2(const Vector& mesh);

/*
 * Function: solve
 * -------------------------------------------
 * Solves a Second order BVP of the form
 * U''(x) = f(x), with the given right hand side
 * and boundary condition, over the given mesh.
 * Receives:
 *  rhs_f   -> Pointer to function that represents
 *              the right hand side.
 *  l_BC    -> constante Left Boundary condition.
 *  r_BC    -> constant Right Boundary Condition.
 *  mesh    -> Vector of uniformly spaced points
 *              between 0 and 1.
 *
 * Saves the solution in a file.
*/
void solve(double (*rhs_f) (double x),
           double l_BC, double r_BC,
           const Vector& mesh);
/*
 * Function: solve
 * -------------------------------------------
 * Solves a Second order BVP of the form
 * U''(x) = f(x), with the given right hand side
 * and boundary condition, over the given mesh.
 * Receives:
 *  rhs_f   -> Pointer to function that represents
 *              the right hand side.
 *  l_BC    -> constante Left Boundary condition.
 *  r_BC    -> constant Right Boundary Condition.
 *  mesh    -> Vector of uniformly spaced points
 *              between 0 and 1.
 *  result  -> Empty vector that will hold the result.
 *  save    -> boolean flag giving the option of
 *              saving the solution on files.
 *
 * Overwrites "result" to contain the approximation, saves the
 * solution on file if "save" flag is true.
*/
void solve(double (*rhs_f) (double x),
           double l_BC, double r_BC,
           const Vector& mesh,
           Vector& result,
           bool save=false);

// Implementation of the case in part e
// of Coursework P1.
void case1();

// Implementation of the case in part d
// of Coursework P1.
void case2();

int main(){

    case1();
    case2();

    return 0;
}

void case1(){

    int n_points=9;
    double alpha = 0;   // Left Boundary Condition.
    double beta = 1;    // Right Boundary Condition.

    // A unifromly spaced mesh between 0 and 1 of n_points+1 points.
    Vector mesh = hlp.linspace(0,1,n_points+1);

    // Right hand side function:
    double (*this_rhs) (double x);
    this_rhs = &rhs_f1;

    // Algorithmic Implementation
    Vector sol(n_points+1);
    solve(this_rhs, alpha, beta, mesh, sol, true);

    Vector true_sol = true_solution1(mesh);
    hlp.saveVectorToFile(true_sol, COMPLETE_PATH+"true.dat");

    // Code Verification
    std::ofstream this_file;
    hlp.openOutputFile(COMPLETE_PATH+"errors_case1.dat",
                       this_file);
    std::cout << "" << std::endl;
    std::cout << "Convergence Analysis for Case1:" << std::endl;
    double mesh_size = 0.5;
    while(mesh_size > 1e-3){
        mesh_size *= 0.5;
        std::cout << "For h= " << mesh_size;
        int this_n = 1/mesh_size;
        Vector this_mesh = hlp.linspace(0,1,this_n);
        Vector this_sol(this_n);
        solve(this_rhs, alpha, beta, this_mesh, this_sol);
        Vector this_true = true_solution1(this_mesh);
        Vector difference = this_true-this_sol;
        double norm = difference.InfinityNorm();
        std::cout << ", error = " << norm << std::endl;
        hlp.saveData(mesh_size, norm, this_file);
    }
    hlp.closeOutputFile(this_file);
}

void case2(){
    double alpha = 0;   // Left Boundary Condition.
    double beta = 1;    // Right Boundary Condition.

    // Right hand side function:
    double (*this_rhs) (double x);
    this_rhs = &rhs_f2;

    // Code Verification
    std::ofstream this_file;
    hlp.openOutputFile(COMPLETE_PATH+"errors_case2.dat",
                       this_file);
    std::cout << "" << std::endl;
    std::cout << "Convergence Analysis for Case2:" << std::endl;
    double mesh_size = 0.5;
    while(mesh_size > 1e-3){
        mesh_size *= 0.5;
        std::cout << "For h= " << mesh_size;
        int this_n = 1/mesh_size;
        Vector this_mesh = hlp.linspace(0,1,this_n);
        Vector this_sol(this_n);
        solve(this_rhs, alpha, beta, this_mesh, this_sol);
        Vector this_true = true_solution2(this_mesh);
        Vector difference = this_true-this_sol;
        double norm = difference.InfinityNorm();
        std::cout << ", error = " << norm << std::endl;
        hlp.saveData(mesh_size, norm, this_file);
    }
    hlp.closeOutputFile(this_file);
}

// Solve function that overwrites result with solution, and saves on file
// depending on "save" flag
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

double true_sol1(double x){
    return x - sin(M_PI*x);
}

double true_sol2(double x){
    return x*sqrt(x);
}

double rhs_f1(double x){
    //return pow(M_PI, 2)*sin(M_PI*x);
    return -50.0/3.0;
}

double rhs_f2(double x){
    return 3.0/(4.0*sqrt(x));
}

Vector true_solution1(const Vector& xvec){
    int n = xvec.GetSize();
    Vector result(n);
    for (int i=0; i<n; i++){
        result[i] = true_sol1(xvec.Read(i));
    }
    return result;
}

Vector true_solution2(const Vector& xvec){
    int n = xvec.GetSize();
    Vector result(n);
    for (int i=0; i<n; i++){
        result[i] = true_sol2(xvec.Read(i));
    }
    return result;
}
