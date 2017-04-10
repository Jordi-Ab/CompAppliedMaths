#ifndef ITERATIVESOLVER_HPP
#define ITERATIVESOLVER_HPP

#include "Matrix.hpp"
#include "Vector.hpp"
#include <stdexcept>

class IterativeSolver{

public:

    // Matrix Constructor
    IterativeSolver(const Matrix& A, const Vector& b,
                    const Vector& initial_x);

    // Tridiagonal Constructor
    IterativeSolver(const Vector& ld,
                    const Vector& d,
                    const Vector& ud,
                    const Vector& b,
                    const Vector& initial_x);

    // Destructor
    ~IterativeSolver();

    // Sets the desired error for the result.
    void setDesiredError(double tol);

    // Sets the tolerated iterations before deciding
    // a failure in convergence.
    void setIterationsTolerance(int tol);

    // Solves the system Ax = b, using Gauss Seidel
    // iterative procedure. Overwrites the given
    // vector with the result.
    void solveGaussSeidel(Vector& result);

    // Solves the system Ax = b, using SOR
    // Succesive Over-Relaxation iterative
    // procedure, with the given omega
    // relaxation parameter. Overwrites the given
    // vector with the result.
    void solveSOR(double omega, Vector& result);

private:

    // Matrix A of the system Ax = b
    Matrix* _A;

    // Vector b (right hand side) of
    // the system Ax = b.
    Vector* _b;

    // Last computed value.
    Vector* _init_guess;

    // Tridiagonal representation of the
    // Matrix A of the system Ax = b,
    // stored in vectors.
    Vector* _ud;
    Vector* _d;
    Vector* _ld;

    //Dimensions of the system.
    int _size;

    // tolerated iterations before deciding
    // a failure in convergence.
    double _iter_tol = 500;

    // desired error for the result.
    double _error_tol = 1e-8;

    // Gauss Seidel step for the complete matrix case
    void gsStep(Vector& x, double omega=-1);

    // Gauss Seidel step for the tridiagonal matrix case
    void gsStepT(Vector& x, double omega=-1);

    // Determines convergence of the method.
    bool hasConverged(Vector& x);

};

/*
 * Function: tridiagTimesVec
 * ------------------------
 * Matrix vector multiplication
 * where the tridiagonal of the matrix
 * is stored on the vectors ld, d, ud
 * for the upper diagonal, main diagonal,
 * and upper diagonal respectiveley.
 */
Vector tridiagTimesVec(const Vector& ld,
                       const Vector& d,
                       const Vector& ud,
                       const Vector& v);

#endif // ITERATIVESOLVER_HPP
