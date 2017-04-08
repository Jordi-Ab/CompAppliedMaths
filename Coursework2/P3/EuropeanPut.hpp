#ifndef EUROPEANPUT_HPP
#define EUROPEANPUT_HPP

#include "AbstractBlackScholes.hpp"

class EuropeanPut:public AbstractBlackScholes{

public:

    /*
     * Constructor: EuropeanPut
     * -------------------------------
     * Sets up everything for AbstractBlackScholes
     * object, and for this object.
     */
    EuropeanPut(int Nx, int Nt, double r,
                double sigma, double T,
                double K, double R);

    /*
     * Function: computeTrueSolution
     * -------------------------------
     * Evaluates the true solution of the European Put
     * at each space point and for all times.
     * Saves the results on an output file named
     * "output_file_name".
     */
    void computeTrueSolution(std::string output_file_name);

    /*
     * Function: solveSemiImplicit
     * -------------------------------
     * computes an approximate solution of the European Put
     * using a Semi Implicit method at each time step, and
     * for all times. Saves the results on an output file named
     * "output_file_name".
     * If compute_error flag is true, it doesn't save the
     * result on a file and calls the computeError function
     * with the value of the solution at the last time step.
     */
    void solveSemiImplicit(std::string output_file_name = "output.dat",
                           bool compute_error=false);

    /*
     * Function: solveFullyImplicit
     * -------------------------------
     * computes an approximate solution of the European Put
     * using a Fully Implicit method at each time step, and
     * for all times. Saves the results on an output file named
     * "output_file_name".
     * If compute_error flag is true, it doesn't save the
     * result on a file and calls the computeError function
     * with the value of the solution at the last time step.
     */
    void solveFullyImplicit(std::string output_file_name= "output.dat",
                            bool compute_error=false);

    /*
     * Function: solveExplicit
     * -------------------------------
     * computes an approximate solution of the European Put
     * using an Explicit method at each time step, and
     * for all times. Saves the results on an output file named
     * "output_file_name".
     * If compute_error flag is true, it doesn't save the
     * result on a file and calls the computeError function
     * with the value of the solution at the last time step.
     */
    void solveExplicit(std::string output_file_name= "output.dat",
                       bool compute_error=false);

    /*
     * Function: getErrorValue
     * -------------------------------
     * Returns the error value of the last solve performed.
     * Error value is saved on the _error_value instance.
     * A solve method with the compute error flag on must
     * first be called.
     */
    double getErrorValue();


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
    double f0(double t);

    // Right Boundary Condition Function fR(t)
    double fR(double t);

    // Initial Condition Function g(x)
    double g(double x);

    // True Solution Function
    double true_sol(double t, double x);

    // Functions needed for the true solution.
    double d1(double t, double x);
    double d2(double t, double x);
    double stdCumNormal(double d); // Cumulative Standard Normal Distribution.

    // Fills in the values for _main_diag, _upp_diag
    // and _low_diag according to a Semi Implicit Discretization.
    void semiImplicitDiscretize(double alpha, double beta);

    // Fills in the values for _main_diag, _upp_diag
    // and _low_diag according to a Fully Implicit Discretization.
    void fullyImplicitDiscretize(double alpha, double beta);

    // Fills in the values for _main_diag, _upp_diag
    // and _low_diag according to an Explicit Discretization.
    void explicitDiscretize(double alpha, double beta);

    /*
     * Function: computeError
     * ---------------------------------------------
     * Computes the Inifinity Norm of the difference
     * between the true solution and the computed solution
     * at the final time step "uh_finalT". This vector is
     * given by the solve function when "compute_error"
     * flag is on.
     */
    void computeError(Vector uh_finalT);


};

#endif // EUROPEANPUT_HPP
