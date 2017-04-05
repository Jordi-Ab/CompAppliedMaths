/************************************************
 * in this script we'll solve the Heat Equation
 * with Dirichlet Boundary Conditions, and one
 * initial condition.
 * nameley:
 *
 * du/dt - a * d^2u/dx^2 = f(x,t)
 *  over: 0 <= x <= 1, and 0 <= t <= T,
 *  with BC's:
 *      u(0, t) = g0
 *      u(1, t) = g1
 *  Initial Condition:
 *      u(x, 0) = u0(x)
 *
 * using an Implicit method.
 ***********************************************/
#include <iostream>
#include <string>
#include <cmath>
#include "Helpers.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "EuropeanPut.hpp"

static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/CompAppliedMaths_git/Coursework2/P3/OutputData/";
Helpers hlp;


void case1();
void case2();
void convergenceAnalysis(EuropeanPut& an_option);

int main(){

    //case1();
    case2();

    return 0;
}

void case1(){
    int Nx = 50; // Number of points in the x dimension.
    int Nt = 500; // Number of time points.
    int T = 5; // Finish time.

    // Parameters of Black Scholes Equation.
    double r = 0; // Interest rate
    double sigma = 0.5; // Volatility
    double K = 100; // Strike Price
    double R = 3*K; // Artificial limit of asset price.

    EuropeanPut option(Nx, Nt, r, sigma, T, K, R);

    std::string U_file_name = COMPLETE_PATH + "true_case1.dat";
    std::string Uh_file_name = COMPLETE_PATH + "computed_case1.dat";
    std::string grid_file_name = COMPLETE_PATH + "grid_case1.dat";

    option.saveGrid(grid_file_name);
    option.computeTrueSolution(U_file_name);
    option.solveSemiImplicit(Uh_file_name);

    // CONVERGENCE ANALYSIS
    convergenceAnalysis(option);
}

void case2(){
    int Nx = 50; // Number of points in the x dimension.
    int Nt = 500; // Number of time points.
    int T = 5; // Finish time.

    // Parameters of Black Scholes Equation.
    double r = 0.1; // Interest rate
    double sigma = 0.1; // Volatility
    double K = 100; // Strike Price
    double R = 3*K; // Artificial limit of asset price.

    EuropeanPut option(Nx, Nt, r, sigma, T, K, R);

    std::string U_file_name = COMPLETE_PATH + "true_case2.dat";
    std::string Uh_file_name = COMPLETE_PATH + "computed_case2.dat";
    std::string grid_file_name = COMPLETE_PATH + "grid_case2.dat";

    option.saveGrid(grid_file_name);
    option.computeTrueSolution(U_file_name);
    option.solveSemiImplicit(Uh_file_name);

    // CONVERGENCE ANALYSIS
    convergenceAnalysis(option);
}

void convergenceAnalysis(EuropeanPut& an_option){
    double C = 4;
    double h = 0.5;
    std::ofstream this_file;
    hlp.openOutputFile(COMPLETE_PATH+"errors.dat", this_file);

    std::cout << "" << std::endl;
    std::cout << "Errors for more refined step sizes h, and time steps delta_t,";
    std::cout << "\ndelta_t is taken to be "<<C<<"*h^2."<< std::endl;
    std::cout << "" << std::endl;
    while(h > 1e-2){
        h *= 0.5;
        double this_dt = C*(h*h);
        int this_Nx = an_option.getAssetPriceLimit()/h;
        int this_Nt = an_option.getFinalTime()/this_dt;

        an_option.setNumSpacePoints(this_Nx);
        an_option.setNumTimePoints(this_Nt);

        an_option.solveSemiImplicit("", true);
        double error = an_option.getErrorValue();

        std::cout << "For h= " << h <<", dt = " << this_dt;
        std::cout << " => error= " << error << std::endl;
        hlp.saveData(h, error, this_file);
    }
    hlp.closeOutputFile(this_file);
}

