#include <iostream>
#include <string>
#include <cmath>
#include "Helpers.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "AmericanPut.hpp"
#include "ProjectedSOR.hpp"

static const std::string COMPLETE_PATH= "/Users/user/Documents/Maestria/Computational Applied Maths/CompAppliedMaths_git/Coursework2/P5/OutputData/";
Helpers hlp; // Set of useful functions

void case1();

int main(){

    case1();
    return 0;
}

void case1(){
    int Nx = 50;            // Number of points in the x dimension.
    int Nt = 500;           // Number of time points.
    int T = 5;              // Finish time.

    // Parameters of Black Scholes Equation.
    double r = 0.05;         // Interest rate
    double sigma = 0.5;     // Volatility
    double K = 100;         // Strike Price
    double R = 3*K;         // Artificial limit of asset price.

    // Instantiate American Put object
    AmericanPut option(Nx, Nt, r, sigma, T, K, R);

    // Output file names
    std::string Uh_file_name = COMPLETE_PATH + "computed_case1.dat";
    std::string grid_file_name = COMPLETE_PATH + "grid_case1.dat";

    // Solve
    option.saveGrid(grid_file_name);
    option.solveSemiImplicit(Uh_file_name);
}
