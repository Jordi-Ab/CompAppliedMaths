#include "ProjectedSOR.hpp"

ProjectedSOR::ProjectedSOR(const Matrix& B, const Vector& f,
                           const Vector& psi, double omega,
                           const Vector& initial_x){

    // Assert Dimesions of given objects are correct.
    if(B.GetNumberOfRows() != B.GetNumberOfCols()){
        throw std::runtime_error("Error on Constructor, not a Square Matrix.");
    }else{
        _size = B.GetNumberOfRows();
    }if(f.GetSize() != _size ||
        psi.GetSize() != _size ||
        initial_x.GetSize() != _size){
        throw std::runtime_error("Error on Constructor. Wrong dimension of vector.");
    }

    _B = new Matrix(B);
    _f = new Vector(f);
    _psi = new Vector(psi);
    _omega = omega;
    _last_x = new Vector(initial_x);

}

ProjectedSOR::~ProjectedSOR(){
    delete _B;
    delete _f;
    delete _psi;
    delete _last_x;
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
    if(result.GetSize() != _size)
        throw std::runtime_error("Dimensions error on solve function.");

    int iterations = 0;
    while (true){
        // Gauss Seidel step. Overwrites result vector.
        gsStep(result);
        iterations += 1;
        Vector difference = result - (*_last_x);
        if (difference.Norm(2) <= _error_tol){
            std::cout << "SOR Method converged in ";
            std::cout << iterations << " iterations." << std::endl;
            break;
        }else if (iterations >= _iter_tol){
            std::cout << "SOR Method failed to converge after ";
            std::cout << iterations << " iterations." << std::endl;
            break;
        }
        // difference between result and last x represents the direction.
        // THE PROJECT POINT IF VIOLATE CONSTRAINT SHOULD GO IN THE GS STEP?
        (*_last_x) = (*_last_x) + difference*_omega; // Perform another iteration.
    }
}

/*
* Criterion to determine convergence of
* a method. I don't use it because involves
* a matrix vector multiplication, vector
* susbtraction and then computing the norm
* of that. */
bool ProjectedSOR::hasConverged(Vector& x){
    Vector residual = (*_B) * x - (*_f);
    double norm = residual.Norm(2);
    if(norm <= _error_tol) return true;
    return false;
}

void ProjectedSOR::gsStep(Vector& x){
    for (int r=0; r<_size; r++){
        x[r] = (*_f)[r];
        for (int c=0; c<r; c++){
            x[r] -= (*_B)(r+1,c+1)*x[c];
        }
        for (int c=r+1; c<_size; c++){
            x[r] -= (*_B)(r+1,c+1) * (*_last_x)[c];
        }
        // Update solution
        //double next_one = x[r] / (*_B)(r+1, r+1);
        x[r] = x[r] / (*_B)(r+1, r+1);
        x[r] = min( ( (*_psi)[r], x[r]);
    }
}
