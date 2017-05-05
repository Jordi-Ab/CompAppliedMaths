#include "AmericanPut.hpp"
#include "ProjectedSOR.hpp"

AmericanPut::AmericanPut(int Nx, int Nt, double r,
                         double sigma, double T,
                         double K, double R){

    // Initialize variables for Black Scholes
    _Nx = Nx;
    _Nt = Nt;
    _r = r;
    _sigma = sigma;
    _T = T;
    _K = K;
    _R = R;

    // Initialize grids on Black Scholes.
    setSpaceGrid();
    setTimeGrid();

}

void AmericanPut::semiImplicitDiscretize(double alpha, double beta){

    _low_diag = new Vector(_Nx-2);
    _main_diag = new Vector(_Nx-1);
    _upp_diag = new Vector(_Nx-2);

    for (int n=1; n<_Nx-1; n++){
        (*_low_diag)[n-1] = -0.5*alpha*( pow(n+1,2.0) ); // Lower diagonal vector
        (*_main_diag)[n-1] = 1 + beta + beta*n + alpha*(n*n); // Main diagonal vector
        (*_upp_diag)[n-1] = -( beta*n + 0.5*alpha*(n*n) ); // Upper diagonal vector.
    }
    (*_main_diag)[-1] = 1 + beta + beta*(_Nx-1) + alpha*( pow(_Nx-1,2.0) );
}

/*
 * Function: solveSemiImplicit
 * -------------------------------
 * computes an approximate solution of the American Put
 * using a Semi Implicit method at each time step, and
 * for all times. Saves the results on an output file named
 * "output_file_name".
 * If compute_error flag is true, it doesn't save the
 * result on a file and calls the computeError function
 * with the value of the solution at the last time step.
 */
void AmericanPut::solveSemiImplicit(std::string output_file_name){
    std::ofstream out_file;

    hlp.openOutputFile(output_file_name, out_file);

    // Constants for discretizing.
    double alpha = pow(_sigma, 2)*_dt;
    double beta = _r*_dt;

    // Outer values of discretized matrix,
    // to be taken from solution:
    double l_0 = -0.5*alpha;/*
    First value of outer lower diagonal*/
    double u_R = -( beta*(_Nx-1) + 0.5*alpha*(pow(_Nx-1,2.0) ));/*
    Last value of outer upper diagonal*/

    // Compute values for inner: lower diagonal, main diagonal,
    // and upper diagonal.
    semiImplicitDiscretize(alpha, beta);

    Vector uh(_Nx-1); // Inner vector of states.
    Vector psi(_Nx-1); // Vector of constraints
    Vector init_x(_Nx-1); // Initial guess
    double omega = 1.8; // Damping factor for Projected SOR.

    // Set first state as Initial Condition,
    // and set vector of constraints at each grid point.
    for (int i=0; i<_Nx-1; i++){
        double x = _s_grid->Read(i+1);
        uh[i] = initCondF(x);
        psi[i] = constraintF(x);
    }

    // Advance forward in time
    for(double t=0; t<=_T; t+=_dt){

        double l_BC = leftBcF(t+_dt); // Left Boundary Condition
        double r_BC = rightBcF(t+_dt); // Right Boundary Condition

        uh[0] = uh[0] - l_0*l_BC;
        uh[-1] = uh[-1] - u_R*r_BC;

        // Instantiate a Projected SOR Solver
        ProjectedSOR psor(*_low_diag, *_main_diag, *_upp_diag,
                          uh, psi, omega, init_x);

        // Solve the system (uh gets overwritten)
        psor.solve(uh);

        // Save solution in file including Boundary Conditions.
        hlp.saveSolution(t, l_BC, uh, r_BC, out_file);
    }

    hlp.closeOutputFile(out_file);
}

// Left Boundary Condition Function f0(t)
double AmericanPut::leftBcF(double t){
    return _K + 0*t;
}

// Right Boundary Condition Function fR(t)
double AmericanPut::rightBcF(double t){
    return 0*t;
}

// Initial Condition Function g(x)
double AmericanPut::initCondF(double x){
    return fmax(_K - x, 0);
}

// Initial Condition Function g(x)
double AmericanPut::constraintF(double x){
    return fmax(_K - x, 0);
}
