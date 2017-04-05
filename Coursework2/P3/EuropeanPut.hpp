#ifndef EUROPEANPUT_HPP
#define EUROPEANPUT_HPP

#include "AbstractBlackScholes.hpp"

class EuropeanPut:public AbstractBlackScholes{
public:

    EuropeanPut(int Nx, int Nt, double r,
                double sigma, double T,
                double K, double R);

    void computeTrueSolution(std::string output_file_name);

    void solveSemiImplicit(std::string output_file_name = "output.dat",
                           bool compute_error=false);

    void solveFullyImplicit(std::string output_file_name);

    void solveExplicit(std::string output_file_name);

    double getErrorValue();


private:

    // Instance Variables
    Vector* _main_diag;
    Vector* _upp_diag;
    Vector* _low_diag;

    double _error_value=-1;

    // Left Boundary Condition Function f0(t)
    double f0(double t);

    // Right Boundary Condition Function fR(t)
    double fR(double t);

    // Initial Condition Function g(x)
    double g(double x);

    // True Solution Function
    double true_sol(double t, double x);

    double d1(double t, double x);
    double d2(double t, double x);
    double stdCumNormal(double d);

    void semiImplicitDiscretize(double alpha, double beta);
    void fullyImplicitDiscretize(double alpha, double beta);
    void explicitDiscretize(double alpha, double beta);

    void computeError(Vector uh_finalt);


};

#endif // EUROPEANPUT_HPP
