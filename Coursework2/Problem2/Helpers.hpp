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
    Helpers();
    double infinityNorm(const Vector& v);
    void printV(const Vector& v);
    bool isAlmostZero(double number, double error_tol=1e-6);
    void saveVectorToFile(const Vector& v, std::string fileName);
    Vector linspace(double first, double last, int n);
    void solveTridiagonal(const Vector& diag, const Vector& offd, Vector& x);
    Vector tridiagTimesVec(const Vector& ld,const Vector& d,
                           const Vector& ud,const Vector& v);
    void openOutputFile(std::string file_name, std::ofstream& output_file);
    void saveData(const double h, double norm, std::ofstream& output_file);
    void saveSolution(const double t, const double l_BC,
                      const Vector& us, const double r_BC,
                      std::ofstream& output_file);
    void saveData(const double t, const Vector& us, std::ofstream& output_file);
    void closeOutputFile(std::ofstream& output_file);

};

#endif // HELPERS_HPP
