#ifndef PROJECTEDSOR_HPP
#define PROJECTEDSOR_HPP

#include "Matrix.hpp"
#include "Vector.hpp"
#include <stdexcept>

class ProjectedSOR{
public:

    // Constructor when given the complete Matrix.
    ProjectedSOR(const Matrix& B, const Vector& f,
                 const Vector& psi, double omega,
                 const Vector& initial_x);

    // Constructor when given the tridiagonal vectors.
    ProjectedSOR(const Vector& low_diag,
                 const Vector& diag,
                 const Vector& upp_diag,
                 const Vector& f,
                 const Vector& psi,
                 double omega,
                 const Vector& initial_x);

    // Destructor
    ~ProjectedSOR();

    // Sets the desired error for the result.
    void setDesiredError(double tol);

    // Sets the tolerated iterations before deciding
    // a failure in convergence.
    void setIterationsTolerance(int tol);

    // Solve the inequality problem
    void solve(Vector& result);

private:

    // Discretized Matrix of the system.
    Matrix* _B;

    // Tridiagonal version of the Discretized
    // Matrix. Stored using vectors
    // for efficiency.
    Vector* _low_diag; // Lower Diagonal
    Vector* _main_diag; // Main Diagonal
    Vector* _upp_diag; // Upper Diagonal

    // Right Hand side vector of the system.
    Vector* _f;

    // Vector of constraints.
    Vector* _psi;

    // Vector of initial guess.
    Vector* _init_guess;

    // Damping factor.
    double _omega;

    //Dimensions of the system.
    int _size;

    // Tolerated iterations before deciding
    // a failure in convergence.
    double _iter_tol = 500;

    // Desired error for the result.
    double _error_tol = 1e-8;

    // Perform one Iteration of SOR Algorithm
    // using complete Matrix _B.
    void makeStep(Vector& result);

    // Perform one Iteration of SOR Algorithm
    // using tridiagonal vectors _low_diag,
    // _main_diag, _upp_diag.
    void makeStepTri(Vector& x);

    // Test for Convergence.
    bool hasConverged(Vector& previous, Vector& next);

};

#endif // PROJECTEDSOR_HPP
