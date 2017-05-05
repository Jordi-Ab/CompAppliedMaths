#include "ProjectedSOR.hpp"
#include "Helpers.hpp"

// Constructor when given the complete Matrix.
ProjectedSOR::ProjectedSOR(const Matrix& B, const Vector& f,
                           const Vector& psi, double omega,
                           const Vector& initial_x){

    // Assert Dimesions of given objects are correct.
    if(B.GetNumberOfRows() != B.GetNumberOfColumns()){
        throw std::runtime_error("Error on Constructor: Not a Square Matrix.");
    }else{
        _size = B.GetNumberOfRows();
    }if(f.GetSize() != _size ||
        psi.GetSize() != _size ||
        initial_x.GetSize() != _size){
        throw std::runtime_error("Error on Constructor: Wrong dimension of vector.");
    }

    _B = new Matrix(B);
    _low_diag = NULL;
    _main_diag = NULL;
    _upp_diag = NULL;

    _f = new Vector(f);
    _psi = new Vector(psi);
    _init_guess = new Vector(initial_x);

    if(omega < 0 || omega > 2){
        throw std::runtime_error("Error on Constructor: Invalid valid for Damping Factor Omega.");
    }
    _omega = omega;


}

// Constructor when given the tridiagonal vectors.
ProjectedSOR::ProjectedSOR(const Vector& low_diag,
                           const Vector& diag,
                           const Vector& upp_diag,
                           const Vector& f,
                           const Vector& psi,
                           double omega,
                           const Vector& initial_x){

    _size = diag.GetSize();
    if(low_diag.GetSize() != _size - 1 ||
        upp_diag.GetSize() != _size - 1 ||
        f.GetSize() != _size ||
        initial_x.GetSize() != _size){
        throw std::runtime_error("Error on Constructor: Wrong dimension of vector.");
    }

    _B = NULL;
    _low_diag = new Vector(low_diag);
    _main_diag = new Vector(diag);
    _upp_diag = new Vector(upp_diag);

    _f = new Vector(f);
    _init_guess = new Vector(initial_x);
    _psi = new Vector(psi);

    if(omega < 0 || omega > 2){
        throw std::runtime_error("Error on Constructor: Invalid valid for Damping Factor Omega.");
    }
    _omega = omega;
}

ProjectedSOR::~ProjectedSOR(){
    delete _B;
    delete _f;
    delete _psi;
    delete _init_guess;
}

// Sets the desired error for the result.
void ProjectedSOR::setDesiredError(double tol){
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
void ProjectedSOR::setIterationsTolerance(int tol){
    if(tol <= 0)
        throw std::runtime_error("Invalid value for iterations tolerance.");
    _iter_tol = tol;
}

void ProjectedSOR::solve(Vector& result){

    // Assert dimension is correct.
    if(result.GetSize() != _size)
        throw std::runtime_error("Dimensions error on solve function.");

    int iterations = 0;

    Vector u0 = (*_init_guess); // Previous
    result = (*_init_guess); // Next one.
    bool converged = false;
    while (true){

        // Gauss Seidel Step with Damping factor.
        // result gets modified.
        if(_B != NULL) makeStep(result); // Whole matrix case
        else makeStepTri(result); // Tridiagonal vectors case.

        iterations += 1;

        // Check convergence
        if (iterations >= _iter_tol){
            // convergence failure.
            result = *_init_guess;
            break;
        }else if(hasConverged(u0, result)){
            converged = true;
            break;
        }else{
            // Make another iteration.
            u0 = result;
        }
    }

    if(converged){
        std::cout << "PSOR Method converged in ";
        std::cout << iterations << " iterations." << std::endl;
    }else{
        std::cout << "PSOR Method failed to converge after ";
        std::cout << iterations << " iterations." << std::endl;
    }
}

bool ProjectedSOR::hasConverged(Vector& previous, Vector&next){
    Vector residual = previous - next;
    double norm = residual.Norm(2);
    if(norm <= _error_tol) return true;
    return false;
}

void ProjectedSOR::makeStep(Vector& x){
    for (int i=1; i<=_size; i++){
        double y = 0.0;
        // Sum values before diagonal
        for (int j=1; j<=i-1; j++){
            y += (*_B)(i,j) * x(j);
        }
        // Sum values after diagonal
        for (int j=i+1; j<=_size; j++){
            y += (*_B)(i,j) * x(j);
        }
        y =  ( (*_f)(i) - y )/(*_B)(i,i);
        // Project point if violates constraint.
        x(i) = std::max( x(i) + _omega*( y - x(i) ), (*_psi)(i) );
    }
}

void ProjectedSOR::makeStepTri(Vector& x){
    for (int i=1; i<=_size; i++){
        double y = 0.0;
        // Sum values before diagonal (lower diagonal)
        if(i != 1){
            y += (*_low_diag)(i-1) * x(i-1);
        }
        // Sum values after diagonal (upper diagonal)
        if(i !=_size){
            y += (*_upp_diag)(i) * x(i+1);
        }
        y =  ( (*_f)(i) - y )/(*_main_diag)(i);
        // Project point if violates constraint.
        x(i) = std::max( x(i) + _omega*( y - x(i) ), (*_psi)(i) );
    }
}
