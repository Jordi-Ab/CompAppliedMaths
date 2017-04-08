#ifndef PROJECTEDSOR_HPP
#define PROJECTEDSOR_HPP

#include "Matrix.hpp"
#include "Vector.hpp"
#include <stdexcept>

class ProjectedSOR
{
public:

    ProjectedSOR(const Matrix& B, const Vector& f,
                 const Vector& psi, double omega,
                 const Vector& initial_x);

    ~ProjectedSOR();

    // Sets the desired error for the result.
    void setDesiredError(double tol);

    // Sets the tolerated iterations before deciding
    // a failure in convergence.
    void setIterationsTolerance(int tol);

    // Solves the system Ax = b, overwrites the given
    // vector with the result.
    void GsSolve(Vector& result);

    void SorSolve(Vector& result);

private:

    Matrix* _B;
    Vector* _f;

    Vector _psi;

    Vector* _last_x;

    double _omega;

    //Dimensions of the system.
    int _size;

    // tolerated iterations before deciding
    // a failure in convergence.
    double _iter_tol = 500;

    // desired error for the result.
    double _error_tol = 1e-8;

    // Gauss Seidel step.
    void gsStep(Vector& result);

    bool hasConverged(Vector& x);

};

#endif // PROJECTEDSOR_HPP
