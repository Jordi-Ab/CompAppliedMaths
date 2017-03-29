/************************************************
 * in this script we'll solve the Heat Equation
 * with Dirichlet Boundary Conditions, and one
 * initial condition.
 * nameley:
 *
 * du/dt - a * d^2u/dx^2 = f(t,x)
 *  over: 0 < x < 1, and 0 < t < T
 *  with BC's:
 *      u(t, 0) = g0
 *      u(t, 1) = g1
 *  Initial Condition:
 *      u(0, x) = u0(x)
 *
 * using an Implicit method.
 ***********************************************/
#include <iostream>
#include <string>
#include "Helpers.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

//static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/Coursework2/Problem2/OutputData/";
static const std::string COMPLETE_PATH= "OutputData/";
Helpers hlp; // Useful functions

Vector evaluateF(double (*f) (double x, double t),
                 const Vector& xvec, double t);

void solve(int t0, int T, double a,
           double l_BC, double r_BC,
           double (*u0_f) (double x, double t),
           double (*rhs_f) (double x, double t),
           double (*true_sol_f) (double x, double t),
           const Vector& xmesh,const Vector& tmesh);

double computeError(int t0, int T, double a,
                    double l_BC, double r_BC,
                    double (*u0_f) (double x, double t),
                    double (*rhs_f) (double x, double t),
                    double (*true_sol_f) (double x, double t),
                    double h);

double true_sol_f(double x, double t){
    return exp(-4*t)*sin(2*M_PI*x) + x;
    //return exp(-4*t)*sin(2*M_PI*x) + exp(-t)*sin(M_PI*x);
}

double u0_f(double x, double t){
    return sin(2*M_PI*x) + x + 0*t;
    //return sin(2*M_PI*x) + sin(M_PI*x);
}

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
    double h = 0.5;
    // Code Verification
    std::ofstream this_file;
    hlp.openOutputFile(COMPLETE_PATH+"errors.dat", this_file);
    while(h > 1e-2){
        h *= 0.5;
        double error = computeError(t0, T, a, g0, g1,_u0_f, _rhs_f, _true_sol_f,h);
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

    // For outputing the results.
    std::string mesh_filename = COMPLETE_PATH + "x_mesh.dat";
    std::string solution_filename = COMPLETE_PATH + "output.dat";
    std::string true_sol_filename = COMPLETE_PATH + "true.dat";

    int Nx = xmesh.GetSize()-1;

    double h = xmesh.Read(1) - xmesh.Read(0);
    double dt = tmesh.Read(1) - tmesh.Read(0);   
    Vector interior_xmesh = hlp.linspace(h,1-h, Nx-1); // mesh of interior points.

    double cfl = a * ( dt/pow(h,2.0) ); // Courant-Friedrichs-Loewy coefficient.

    hlp.saveVectorToFile(xmesh, mesh_filename);

    // Sparse Inner Tridiagonal Matrix.
    // Its the sparse representation of K.
    Vector diag(Nx-1);
    Vector offd(Nx-2);
    diag.fillWith(2);
    offd.fillWith(-1);

    Vector eye(Nx-1);
    eye.fillWith(1);

    // ( I + cfl*K ):
    diag = eye + diag*cfl;
    offd = offd*cfl;

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
        Vector us gets modified to contain the new solution.*/
    }
    hlp.closeOutputFile(output_file);
    hlp.closeOutputFile(true_sol_output);
}

double computeError(int t0, int T, double a,
                    double l_BC, double r_BC,
                    double (*u0_f) (double x, double t),
                    double (*rhs_f) (double x, double t),
                    double (*true_sol_f) (double x, double t),
                    double h){

    double dt = 4*(h*h);
    int Nx = 1/h;
    std::cout << "For h= " << h;
    std::cout << ", and dt= " << dt << ": ";

    Vector interior_xmesh = hlp.linspace(h,1-h, Nx-1); // mesh of interior points.

    double cfl = a * ( dt/pow(h,2.0) ); // Courant-Friedrichs-Loewy coefficient.

    // Sparse Inner Tridiagonal Matrix.
    // Its the sparse representation of K.
    Vector diag(Nx-1);
    Vector offd(Nx-2);
    diag.fillWith(2);
    offd.fillWith(-1);

    Vector eye(Nx-1);
    eye.fillWith(1);

    // ( I + cfl*K ):
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
        Vector us gets modified to contain the new solution.*/
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
