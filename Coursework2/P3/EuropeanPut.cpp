#include "EuropeanPut.hpp"

EuropeanPut::EuropeanPut(int Nx, int Nt, double r,
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

void EuropeanPut::solveExplicit(std::string output_file_name,
                                bool compute_error){
    std::ofstream out_file;
    hlp.openOutputFile(output_file_name, out_file);

    // Constants for discretizing.
    double alpha = pow(_sigma, 2)*_dt;
    double beta = _r*_dt;

    explicitDiscretize(alpha, beta);
    for(double t=0; t<=_T; t+=_dt){
        double l_BC = f0(t+_dt);
        double r_BC = fR(t+_dt);


    }

    hlp.closeOutputFile(out_file);
}

void EuropeanPut::solveFullyImplicit(std::string output_file_name,
                                     bool compute_error){
    std::ofstream out_file;
    hlp.openOutputFile(output_file_name, out_file);

    // Constants for discretizing.
    double alpha = pow(_sigma, 2)*_dt;
    double beta = _r*_dt;

    fullyImplicitDiscretize(alpha, beta);
    for(double t=0; t<=_T; t+=_dt){
        double l_BC = f0(t+_dt);
        double r_BC = fR(t+_dt);


    }

    hlp.closeOutputFile(out_file);
}

void EuropeanPut::generalSolve(std::string output_file_name,
                               bool compute_error){
    std::ofstream out_file;
    bool save = true;
    if (compute_error) save=false;

    if (save)
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

    // Set first state as Initial Condition
    for (int i=0; i<_Nx-1; i++){
        double x = _s_grid->Read(i+1);
        uh[i] = g(x);
    }

    // Advance forward in time
    for(double t=0; t<=_T; t+=_dt){

        double l_BC = f0(t+_dt);
        double r_BC = fR(t+_dt);

        uh[0] = uh[0] - l_0*l_BC;
        uh[-1] = uh[-1] - u_R*r_BC;

        // Solve the system (rhs gets overwritten)
        hlp.solveTridiagonal(*_low_diag, *_main_diag, *_upp_diag, uh);

        // Save solution in file including Boundary Conditions.
        if(save) hlp.saveSolution(t, l_BC, uh, r_BC, out_file);
    }

    if(compute_error) computeError(uh);/*
    Compute error with uh at the final time T.*/

    if(save) hlp.closeOutputFile(out_file);


}

void EuropeanPut::solveSemiImplicit(std::string output_file_name,
                                    bool compute_error){
    std::ofstream out_file;
    bool save = true;
    if (compute_error) save=false;

    if (save)
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

    // Set first state as Initial Condition
    for (int i=0; i<_Nx-1; i++){
        double x = _s_grid->Read(i+1);
        uh[i] = g(x);
    }

    // Advance forward in time
    for(double t=0; t<=_T; t+=_dt){

        double l_BC = f0(t+_dt);
        double r_BC = fR(t+_dt);

        uh[0] = uh[0] - l_0*l_BC;
        uh[-1] = uh[-1] - u_R*r_BC;

        // Solve the system (rhs gets overwritten)
        hlp.solveTridiagonal(*_low_diag, *_main_diag, *_upp_diag, uh);

        // Save solution in file including Boundary Conditions.
        if(save) hlp.saveSolution(t, l_BC, uh, r_BC, out_file);
    }

    if(compute_error) computeError(uh);/*
    Compute error with uh at the final time T.*/

    if(save) hlp.closeOutputFile(out_file);
}

void EuropeanPut::computeError(Vector uh_finalT){

    // Compute the true solution
    // at the final time T.
    Vector true_s(_Nx-1);
    for (int i=0; i<_Nx-1; i++){
        double x = _s_grid->Read(i+1);
        true_s[i] = true_sol(_T, x);
    }

    Vector difference = uh_finalT-true_s;
    double norm = difference.InfinityNorm();
    _error_value = norm;
}

void EuropeanPut::semiImplicitDiscretize(double alpha, double beta){

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

void EuropeanPut::fullyImplicitDiscretize(double alpha, double beta){

    _low_diag = new Vector(_Nx-2);
    _main_diag = new Vector(_Nx-1);
    _upp_diag = new Vector(_Nx-2);

    for (int n=1; n<_Nx-1; n++){
        (*_low_diag)[n-1] = 0.5*( beta*(n+1) - alpha*( pow(n+1,2.0) ) ); // Lower diagonal vector
        (*_main_diag)[n-1] = 1 + beta + alpha*(n*n); // Main diagonal vector
        (*_upp_diag)[n-1] = -0.5*( beta*n + alpha*(n*n) ); // Upper diagonal vector.
    }
    (*_main_diag)[-1] = 1 + beta + alpha*( pow(_Nx-1,2.0) );
}

void EuropeanPut::explicitDiscretize(double alpha, double beta){
    throw std::runtime_error(" solveExplicit not yet Implemented ");
}

void EuropeanPut::computeTrueSolution(std::string output_file_name){
    std::ofstream out_file;
    hlp.openOutputFile(output_file_name, out_file);
    for(double t=0; t<=_T; t+=_dt){
        Vector true_s(_Nx+1);
        for (int i=0; i<=_Nx; i++){
            double x = _s_grid->Read(i);
            true_s[i] = true_sol(t, x);
        }
        hlp.saveData(t, true_s, out_file);
    }
    hlp.closeOutputFile(out_file);
}

// Left Boundary Condition Function f0(t)
double EuropeanPut::f0(double t){
    return _K*exp(-_r*t);
    //return _K + 0*t;
}

// Right Boundary Condition Function fR(t)
double EuropeanPut::fR(double t){
    return true_sol(t, _R);
    //return 0*t;
}

// Initial Condition Function g(x)
double EuropeanPut::g(double x){
    return fmax(_K - x, 0);
}

// True Solution Function
double EuropeanPut::true_sol(double t, double x){
    if(t==0) return g(x);
    return _K*exp(-_r*t)*stdCumNormal(-d1(t,x)) - x*stdCumNormal(-d2(t,x));
}

double EuropeanPut::d1(double t, double x){
    return ( log(x/_K) + (_r-0.5*(_sigma*_sigma))*t ) / (_sigma*sqrt(t));
}

double EuropeanPut::d2(double t, double x){
    return ( log(x/_K) + (_r+0.5*(_sigma*_sigma))*t ) / (_sigma*sqrt(t));
}

double EuropeanPut::stdCumNormal(double d){
    return 0.5*std::erfc( -(d/sqrt(2)) );
}

double EuropeanPut::getErrorValue(){
    if (_error_value == -1){
        std::cout << "No error value computed yet.";
        std::cout << "Run a solve function with compute error flag = true";
        std::cout << std::endl;
        return -1;
    }else{
        return _error_value;
    }
}

