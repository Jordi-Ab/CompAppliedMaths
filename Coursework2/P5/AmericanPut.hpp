#ifndef AMERICANPUT_HPP
#define AMERICANPUT_HPP

#include "AbstractBlackScholes.hpp"
#include "AbstractIterativeSolver.hpp"

class AmericanPut:public AbstractBlackScholes{
public:

    /*
     * Constructor: AmericanPut
     * -------------------------------
     * Sets up everything for AbstractBlackScholes
     * object, and for this object.
     */
    AmericanPut(int Nx, int Nt, double r,
                double sigma, double T,
                double K, double R);

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
    void solveSemiImplicit(std::string output_file_name = "output.dat");

private:

    // INSTANCE VARIABLES

    // Main diagonal of discretized tridiagonal matrix.
    Vector* _main_diag;
    // Upper diagonal of discretized tridiagonal matrix.
    Vector* _upp_diag;
    // Lower diagonal of discretized tridiagonal matrix.
    Vector* _low_diag;

    // Will hold the error value of the last solve
    // that had the compute_error flag on.
    double _error_value=-1;

    // Left Boundary Condition Function f0(t)
    double leftBcF(double t);

    // Right Boundary Condition Function fR(t)
    double rightBcF(double t);

    // Initial Condition Function g(x)
    double initCondF(double x);

    // Function of constraints.
    double constraintF(double x);

    // Fills in the values for _main_diag, _upp_diag
    // and _low_diag according to a Semi Implicit Discretization.
    void semiImplicitDiscretize(double alpha, double beta);

    // Default function for the
    // Left Boundary Condition f0(t)
    double defaultLeftBC(double t){
        return _K + 0*t;
    }

    // Default function for the
    // right Boundary Condition fR(t)
    double defaultRightBC(double t){
        return 0*t;
    }

    // Default function for the
    // Initial Condition g(x)
    double defaultInitCond(double x){
        return fmax(_K - x, 0);
    }

};

#endif // AMERICANPUT_HPP
