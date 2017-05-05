#include <iostream>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "ProjectedSOR.hpp"
#include "Helpers.hpp"

static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/CompAppliedMaths_git/Coursework2/P4/OutputData/";
Helpers hlp;

/*
 * Function: psi_f
 * -------------------
 * Function of constraints.
*/
double psi_f(double x);

/*
 * Function: rhs_f
 * -------------------
 * Right hand side function of the system
*/
double rhs_f(double x);

/*
 * Function: vectorEvaluateF
 * -------------------
 * Evaluates the given function "f" at each
 * point inside the vector "xs", stores the
 * result on the vector "result".
 * # Equivalent to a vector evaluation of
 *   a function.
*/
void vectorEvaluateF(double (*f) (double x),
                     const Vector & xs,
                     Vector& result);

// Solves the Elliptic Inequality problem using the
// Projected SOR method.
void solveEllipticInequality(double (*rhs_f) (double x),
                             double (*psi_f) (double x),
                             double l_BC, double r_BC,
                             const Vector& mesh,
                             Vector& result, bool save);

// unconstrained true solution of
// the standard elliptic model
Vector unconstrainedTrueSol(const Vector& xvec);

int main(){

    int Nx=16;
    double l_BC = 0;   // Left Boundary Condition.
    double r_BC = 0;    // Right Boundary Condition.

    // A unifromly spaced mesh between 0 and 1 of Nx+1 points.
    Vector mesh = hlp.linspace(0,1,Nx+1);

    // Right hand side function:
    double (*this_rhs) (double x);
    this_rhs = &rhs_f;

    // Function of constraints.
    double (*this_psi) (double x);
    this_psi = &psi_f;

    // Algorithmic Implementation
    Vector solution(Nx+1);
    solveEllipticInequality(this_rhs, this_psi,
                            l_BC, r_BC, mesh,
                            solution, true);

    Vector true_sol = unconstrainedTrueSol(mesh);
    hlp.saveVectorToFile(true_sol, COMPLETE_PATH+"unconstrained_true.dat");

    return 0;
}

// Solve function that overwrites result with solution, and saves on file
// depending on "save" flag
void solveEllipticInequality(double (*rhs_f) (double x),
                             double (*psi_f) (double x),
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
    Vector psi(inner_n); // Vector of constraints
    Vector init_x(inner_n); // Initial guess

    // Left Boundary condition on Right hand side
    rhs[0] = rhs_f( mesh.Read(1) ) - l_BC/(h*h);

    // Assemble inner tridiagonal "matrix"
    // and inner right hand side.
    for (int i = 0; i<inner_n-1; i++){
        diagonal[i] = -2.0/(h*h);
        off_diag[i] = 1.0/(h*h);
        rhs[i+1] = rhs_f( mesh.Read(i+2) );
        psi[i] = psi_f( mesh.Read(i+1) );
    }
    // Right Boundary condition on Right hand side
    rhs[inner_n-1] = rhs_f( mesh.Read(n-2) ) - r_BC/(h*h);

    diagonal[inner_n-1] = -2.0/(h*h);
    psi[inner_n-1] = psi_f( mesh.Read(n-2) );

    double omega = 1.8; // Damping Factor

    // Instantiate Projected SOR Solver.
    ProjectedSOR solver(off_diag, diagonal, off_diag, rhs, psi, omega, init_x);

    solver.solvePSOR(rhs, COMPLETE_PATH+"iterations.dat");

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

double psi_f(double x){
    return 1 + x*0;
}

double rhs_f(double x){
    return -50.0/3.0 + x*0;
}

void vectorEvaluateF(double (*f) (double x),
                     const Vector & xs,
                     Vector& result){

    for (int i=0; i<xs.GetSize(); i++){
        result[i] = f(xs.Read(i));
    }
}

double unconstrainedTrueSol(double x){
    return (25.0/3.0)*x - (25.0/3.0)*(x*x);
}

Vector unconstrainedTrueSol(const Vector& xvec){
    int n = xvec.GetSize();
    Vector result(n);
    for (int i=0; i<n; i++){
        result[i] = unconstrainedTrueSol(xvec.Read(i));
    }
    return result;
}

