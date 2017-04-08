/************************************************
 * in this script we'll solve the Heat Equation
 * with Dirichlet Boundary Conditions, and one
 * initial condition.
 * nameley:
 *
 * du/dt - a * d^2u/dx^2 = f(x,t)
 *  over: 0 <= x <= 1, and 0 <= t <= T,
 *  with BC's:
 *      u(0, t) = g0
 *      u(1, t) = g1
 *  Initial Condition:
 *      u(x, 0) = u0(x)
 *
 * using an Implicit method.
 ***********************************************/
#include <iostream>
#include <string>
#include "Helpers.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/CompAppliedMaths_git/Coursework2/P2/OutputData/";
//static const std::string COMPLETE_PATH= "OutputData/";
Helpers hlp; // Useful functions

/*
 * Function: evaluateF
 * -------------------------
 * f(x, t) where x is a vector and t is a point.
 * Vector evaluation of the given function "f"
 * at the given Vector "xvec", and double t.
 *
 * Receives:
 *  f -> Pointer to the function that will
 *       get evaluated.
 *  xvec -> Vector of x points where f will
 *       get evaluated.
 *  t -> double, time point.
 *
 * Returns
 *  A Vector containing the results of f(xi, t),
 *  for xi components on xvec.
 */
Vector evaluateF(double (*f) (double x, double t),
                 const Vector& xvec, double t);

/*
 * Function: solve
 * -------------------------
 * Solves a heat equation PDE as the one
 * described at the comments on top,
 * using an implicit method. Saves the
 * solution at each time step on an output file
 *
 * Receives:
 *  t0 -> Starting time point
 *  T -> Final time point
 *  a -> Heat-conduction coefficient of
 *      the Heat equation.
 *  l_BC -> Constant left boundary condition
 *  r_BC -> Constant right boundary condition
 *  u0_f -> Pointer to Initial Condition Function
 *  rhs_f -> Pointer to Right Hand Side Function
 *  true_sol_f -> Pointer to True Solution Function
 *  xmesh -> Equally spaced points between 0 and 1,
 *      for the space dimension.
 *  tmesh -> Equally spaced points between t0 and T,
 *      for the time steps.
 */
void solve(int t0, int T, double a,
           double l_BC, double r_BC,
           double (*u0_f) (double x, double t),
           double (*rhs_f) (double x, double t),
           double (*true_sol_f) (double x, double t),
           const Vector& xmesh,const Vector& tmesh);

/*
 * Function: computeError
 * -------------------------
 * Computes the error between the approximated solution
 * and the true solution at the final time step T, as the
 * infinity norm of the difference U - Uh.
 *
 * Performs all the steps as the above function solve to get
 * the computed solution at the final time T, and then
 * computes the error using Uh @ T.
 *
 * Receives:
 *  # See the above comments for function solve.
 *  h -> Step Size to be used for computing the error.
 *  C -> Constant that couples the step size h with the
 *      time step delta t. i.e. dt = C*h^2
 */
double computeError(int t0, int T, double a,
                    double l_BC, double r_BC,
                    double (*u0_f) (double x, double t),
                    double (*rhs_f) (double x, double t),
                    double (*true_sol_f) (double x, double t),
                    double h, double C);

// True Solution function of the Heat Equation
// intended to be approximated.
double true_sol_f(double x, double t){
    return exp(-4*t)*sin(2*M_PI*x) + x;
    //return exp(-4*t)*sin(2*M_PI*x) + exp(-t)*sin(M_PI*x);
}

// Initial Condition function of the Heat Equation
// intended to be approximated.
double u0_f(double x, double t){
    return sin(2*M_PI*x) + x + 0*t;
    //return sin(2*M_PI*x) + sin(M_PI*x);
}

// Right Hand Side function of the Heat Equation
// intended to be approximated.
double rhs_f(double x, double t){
    return 0*(x+t); // just return a zero.
}

int main(){

    int Nx = 16; // Number of points in the x dimension.
    int Nt = 64; // Number of time points.
    int t0 = 0; // Starting time
    int T = 1; // Finish time.
    double a = 1.0 / pow(M_PI,2.0); // Heat-conduction coefficient.
    double g0 = 0; // Left Boundary Condition.
    double g1 = 1; // Right Boundary Condition.
    double (*_rhs_f) (double x, double t); // Right hand Side function f(x, t)
    _rhs_f = &rhs_f;
    double (*_u0_f) (double x, double t); // Initial Condition Function
    _u0_f = &u0_f;
    double (*_true_sol_f) (double x, double t); // True solution function.
    _true_sol_f = &true_sol_f;

    // Uniformly spaced grid between 0 and 1, for the x variable.
    Vector xmesh = hlp.linspace(0,1, Nx+1);
    // Uniformly spaced grid between t0 and T, for the time variable.
    Vector tmesh = hlp.linspace(t0,T, Nt+1);

    solve(t0, T, a, g0, g1, _u0_f, _rhs_f, _true_sol_f, xmesh, tmesh);

    // CONVERGENCE ANALYSIS
    double C = 4;
    double h = 0.5;
    std::ofstream this_file;
    hlp.openOutputFile(COMPLETE_PATH+"errors.dat", this_file);

    std::cout << "" << std::endl;
    std::cout << "Errors for more refined step sizes h, and time steps delta_t,";
    std::cout << "\ndelta_t is taken to be "<<C<<"*h^2."<< std::endl;
    std::cout << "" << std::endl;
    while(h > 1e-2){
        h *= 0.5;
        double error = computeError(t0, T, a, g0, g1,_u0_f, _rhs_f, _true_sol_f,h, C);
        std::cout << ", error= " << error << std::endl;
        hlp.saveData(h, error, this_file);
    }
    hlp.closeOutputFile(this_file);
    return 0;
}

void solve(int t0, int T, double a,
           double l_BC, double r_BC,
           double (*u0_f) (double x, double t),
           double (*rhs_f) (double x, double t),
           double (*true_sol_f) (double x, double t),
           const Vector& xmesh,const Vector& tmesh){

    // Name of output files.
    std::string mesh_filename = COMPLETE_PATH + "x_mesh.dat";
    std::string solution_filename = COMPLETE_PATH + "output.dat";
    std::string true_sol_filename = COMPLETE_PATH + "true.dat";

    int Nx = xmesh.GetSize()-1;

    double h = xmesh.Read(1) - xmesh.Read(0);
    double dt = tmesh.Read(1) - tmesh.Read(0);   
    Vector interior_xmesh = hlp.linspace(h,1-h, Nx-1); // mesh of interior points.

    double cfl = a * ( dt/pow(h,2.0) ); // Courant-Friedrichs-Loewy coefficient.

    hlp.saveVectorToFile(xmesh, mesh_filename); // save mesh on file.

    // Sparse Inner Tridiagonal Matrix.
    // (Its the sparse representation of K).
    Vector diag(Nx-1);
    Vector offd(Nx-2);
    diag.fillWith(2); // Set all values of vector as 2.
    offd.fillWith(-1); // Set all values of vector as -1.

    // Sparse identity matrix.
    Vector eye(Nx-1);
    eye.fillWith(1); // Set all values of vector as 1.

    // Discretize ( I + cfl*K ):
    diag = eye + diag*cfl;
    offd = offd*cfl;

    // Open files for saving.
    std::ofstream output_file;
    std::ofstream true_sol_output;
    hlp.openOutputFile(solution_filename, output_file);
    hlp.openOutputFile(true_sol_filename, true_sol_output);

    // Advance it forward in time:
    Vector uh = evaluateF(u0_f, interior_xmesh, 0); // Initial Condition
    for(double t=t0; t<=T; t+=dt){

        // Save solution
        hlp.saveSolution(t, l_BC, uh, r_BC, output_file);
        // Compute true solution
        Vector true_s = evaluateF(true_sol_f, xmesh, t);
        // Save true solution
        hlp.saveData(t, true_s, true_sol_output);

        // Assemble Right Hand Side
        Vector rhs = evaluateF(rhs_f, interior_xmesh, t+dt); /*
        Evaluate right hand side function at the interior mesh, and next time.*/
        rhs[0] = rhs[0] + ( a/(h*h) )*l_BC; // Include Left BC at first entry
        rhs[Nx - 2] = rhs[Nx - 2] + ( a/(h*h) )*r_BC; // Include Right BC at last entry

        uh = uh + rhs*dt; // Complete right hand side.

        // Solve the system
        hlp.solveTridiagonal(diag, offd, uh);/*
        Vector uh gets modified to contain the new solution.*/
    }
    hlp.closeOutputFile(output_file);
    hlp.closeOutputFile(true_sol_output);
}

double computeError(int t0, int T, double a,
                    double l_BC, double r_BC,
                    double (*u0_f) (double x, double t),
                    double (*rhs_f) (double x, double t),
                    double (*true_sol_f) (double x, double t),
                    double h, double C){

    double dt = C*(h*h);
    int Nx = 1/h;
    std::cout << "For h= " << h;
    std::cout << ", and dt= " << dt << ": ";

    Vector interior_xmesh = hlp.linspace(h,1-h, Nx-1); // mesh of interior points.

    double cfl = a * ( dt/pow(h,2.0) ); // Courant-Friedrichs-Loewy coefficient.

    // Sparse Inner Tridiagonal Matrix.
    // (Its the sparse representation of K).
    Vector diag(Nx-1);
    Vector offd(Nx-2);
    diag.fillWith(2);
    offd.fillWith(-1);

    Vector eye(Nx-1); // Sparse Identity matrix
    eye.fillWith(1);

    // Discretize ( I + cfl*K ):
    diag = eye + diag*cfl;
    offd = offd*cfl;

    // Advance it forward in time:
    Vector uh = evaluateF(u0_f, interior_xmesh, 0); // Initial Condition
    for(double t=t0; t<=T; t+=dt){
        // Assemble Right Hand Side
        Vector rhs = evaluateF(rhs_f, interior_xmesh, t+dt); /*
        Evaluate right hand side function at the interior mesh, and next time.*/
        rhs[0] = rhs[0] + ( a/(h*h) )*l_BC; // Include Left BC at first entry
        rhs[Nx-2] = rhs[Nx-2] + ( a/(h*h) )*r_BC; // Include Right BC at last entry

        uh = uh + rhs*dt; // Complete right hand side.

        // Solve the system
        hlp.solveTridiagonal(diag, offd, uh);/*
        Vector uh gets modified to contain the new solution.*/
    }

    // Compute the error at the final time T.
    Vector true_s = evaluateF(true_sol_f, interior_xmesh, T);
    Vector difference = uh-true_s;
    double norm = difference.InfinityNorm();
    return norm;
}

// Evaluate a given function at vector x and time t.
Vector evaluateF(double (*f) (double x, double t),
                 const Vector& xvec, double t){
    int n = xvec.GetSize();
    Vector result(n);
    for (int i=0; i<n; i++){
        result[i] = f(xvec.Read(i), t);
    }
    return result;
}
