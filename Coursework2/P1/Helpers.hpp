#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "Vector.hpp"

class Helpers{

public:

    // Default Constructor
    Helpers();

    // Return the infinity norm of a vector.
    double infinityNorm(const Vector& v);

    // Prints a vector on screen in the format:
    // [1, 2, 3, 4]
    void printV(const Vector& v);

    // Test if a number is zero at the given error tolerance.
    bool isAlmostZero(double number, double error_tol=1e-6);

    // Saves the elements of a vector on the given file.
    void saveVectorToFile(const Vector& v, std::string fileName);

    // Vector of n uniformly spaced points between "first" and "last"
    Vector linspace(double first, double last, int n);

    // Solves a tridiagonal system where the upper diagonal and lower diagonal are the same.
    void solveTridiagonal(const Vector& diag, const Vector& offd, Vector& x);

    // Opens an output file and formats it to save data.
    void openOutputFile(std::string file_name, std::ofstream& output_file);

    // Saves data on the output file.
    void saveData(const double h, double norm, std::ofstream& output_file);

    // Closes the given file.
    void closeOutputFile(std::ofstream& output_file);

};

#endif // HELPERS_HPP
