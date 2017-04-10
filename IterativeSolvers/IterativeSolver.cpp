#include "IterativeSolver.hpp"

IterativeSolver::IterativeSolver(const Matrix& A, const Vector& b,
                                 const Vector& initial_x){

    // Assert Dimesions of given objects are correct.
    if(A.GetNumberOfRows() != A.GetNumberOfColumns()){
        throw std::runtime_error("Error on Constructor: Not a Square Matrix.");
    }else{
        _size = A.GetNumberOfRows();
    }if(b.GetSize() != _size ||initial_x.GetSize() != _size){
        throw std::runtime_error("Error on Constructor: Wrong dimension of vector.");
    }

    _A = new Matrix(A);
    _ld = NULL;
    _d = NULL;
    _ud = NULL;
    _b = new Vector(b);
    _init_guess = new Vector(initial_x);
}

IterativeSolver::IterativeSolver(const Vector& ld,
                                 const Vector& d,
                                 const Vector& ud,
                                 const Vector& b,
                                 const Vector& initial_x){

    _size = d.GetSize();
    if(ld.GetSize() != _size - 1 ||
        ud.GetSize() != _size - 1 ||
        b.GetSize() != _size ||
        initial_x.GetSize() != _size){
        throw std::runtime_error("Error on Constructor: Wrong dimension of vector.");
    }

    _A = NULL;
    _ld = new Vector(ld);
    _d = new Vector(d);
    _ud = new Vector(ud);
    _b = new Vector(b);
    _init_guess = new Vector(initial_x);
}

IterativeSolver::~IterativeSolver(){
    delete _A;
    delete _b;
    delete _init_guess;
}

// Sets the desired error for the result.
void IterativeSolver::setDesiredError(double tol){
    if(tol <= 0){
        throw std::runtime_error("Invalid value for desired error.");
    }else if (tol > 1e-10){
        std::cout << "Warning, desired error is too small,";
        std::cout << " might take longer to converge." << std::endl;
    }
    _error_tol = tol;
}

// Sets the tolerated iterations before deciding
// a failure in convergence.
void IterativeSolver::setIterationsTolerance(int tol){
    if(tol <= 0)
        throw std::runtime_error("Invalid value for iterations tolerance.");
    _iter_tol = tol;
}

void IterativeSolver::solveGaussSeidel(Vector& result){
    // Assert dimension is correct.
    if(result.GetSize() != _size)
        throw std::runtime_error("Dimensions error on solve function.");

    int iterations = 0;

    // Initialise solution to initial guess
    result = (*_init_guess);
    bool converged = hasConverged(result);
    while (!converged){
        // Gauss Seidel step.
        if(_A == NULL) gsStepT(result);
        else gsStep(result);
        iterations += 1;
        // Check iterations tolerance.
        if (iterations >= _iter_tol) break;
        // Test Convergence.
        converged = hasConverged(result);
    }

    if(converged){
        std::cout << "Gauss Seidel converged in ";
        std::cout << iterations << " iterations." << std::endl;
    }else{
        std::cout << "Gauss Seidel failed to converge after ";
        std::cout << iterations << " iterations." << std::endl;
    }
}

void IterativeSolver::solveSOR(double omega, Vector& result){
    // Assert dimension is correct.
    if(result.GetSize() != _size)
        throw std::runtime_error("Dimensions error on solve function.");

    if(omega < 0 || omega > 2){
        throw std::runtime_error("Invalid valid for Damping Factor Omega.");
    }

    int iterations = 0;

    // Initialise solution to initial guess
    result = (*_init_guess);
    bool converged = hasConverged(result);
    while (!converged){
        // Gauss Seidel step with Damping factor.
        if(_A == NULL) gsStepT(result, omega);
        else gsStep(result, omega);
        iterations += 1;
        // Check iterations tolerance.
        if (iterations >= _iter_tol) break;
        // Test Convergence.
        converged = hasConverged(result);
    }

    if(converged){
        std::cout << "SOR Method converged in ";
        std::cout << iterations << " iterations." << std::endl;
    }else{
        std::cout << "SOR Method failed to converge after ";
        std::cout << iterations << " iterations." << std::endl;
    }
}

void IterativeSolver::gsStepT(Vector& x, double omega){
    for (int i=1; i<=_size; i++){
        double y = 0.0;
        // Sum values before diagonal
        if(i != 1){
            y += (*_ld)(i-1) * x(i-1);
        }
        // Sum values after diagonal
        if(i !=_size){
            y += (*_ud)(i) * x(i+1);
        }
        y =  ( (*_b)(i) - y )/(*_d)(i);
        if(omega != -1){
            // SOR Step.
            x(i) += omega*( y - x(i) );
        }else{
            // Gauss Seidel step.
            x(i) = y;
        }
    }
}

void IterativeSolver::gsStep(Vector& x, double omega){
    for (int i=1; i<=_size; i++){
        double y = 0.0;
        // Sum values before diagonal
        for (int j=1; j<=i-1; j++){
            y += (*_A)(i,j) * x(j);
        }
        // Sum values after diagonal
        for (int j=i+1; j<=_size; j++){
            y += (*_A)(i,j) * x(j);
        }
        y =  ( (*_b)(i) - y )/(*_A)(i,i);
        if(omega != -1){
            // SOR Step.
            x(i) += omega*( y - x(i) );
        }else{
            // Gauss Seidel step.
            x(i) = y;
        }
    }
}

/*
 * Criterion to determine convergence of
 * a method. I don't use it because involves
 * a matrix vector multiplication, vector
 * susbtraction and then computing the norm
 * of that.
*/
bool IterativeSolver::hasConverged(Vector& x){
    Vector residual(x.GetSize());
    // Is Ax = b yet?
    if(_A != NULL){
        residual = (*_A) * x - (*_b); // Ax - b
    }else{
        Vector mul = tridiagTimesVec(*_ld, *_d, *_ud, x); // Ax
        residual = mul - (*_b); // Ax - b
    }
    double norm = residual.Norm(2);
    if(norm <= _error_tol) return true;
    return false;
}

Vector tridiagTimesVec(const Vector& ld,
                       const Vector& d,
                       const Vector& ud,
                       const Vector& v){
    int n = d.GetSize();
    Vector result(n);
    result[0] = d.Read(0)*v.Read(0) + ud.Read(0)*v.Read(1);
    for (int i=0; i<ld.GetSize()-1; i++){
        result[i+1] = ld.Read(i)*v.Read(i) + d.Read(i+1)*v.Read(i+1) + ud.Read(i)*v.Read(i+2);
    }
    result[n-1] = ld.Read(n-2)*v.Read(n-2) + d.Read(n-1)*v.Read(n-1);
    return result;
}
