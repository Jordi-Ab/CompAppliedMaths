#ifndef ABSTRACTBLACKSCHOLES_HPP
#define ABSTRACTBLACKSCHOLES_HPP

#include <cmath>
#include "Vector.hpp"
#include "Helpers.hpp"
#include <stdexcept>

class AbstractBlackScholes{
public:

    // Getters
    double getNumSpacePoints();
    double getNumTimePoints();
    double getInterestRate();
    double getVolatility();
    double getFinalTime();
    double getStrikePrice();
    double getAssetPriceLimit();
    Vector getSpaceGrid();
    Vector getTimeGrid();
    double getDeltaS();
    double getDeltaT();

    // Setters
    void setNumSpacePoints(int Nx);
    void setNumTimePoints(int Nt);
    void setInterestRate(double r);
    void setVolatility(double sigma);
    void setFinalTime(double T);
    void setStrikePrice(double K);
    void setAssetPriceLimit(double R);
    void setOutputFileName(std::string a_name);
    void setCompletePath(std::string a_path);

    void saveGrid(std::string file_name);



protected:

    Helpers hlp; // Set of useful functions.

    // Instance Variables
    int _Nx;        // Number of points in the space grid.
    int _Nt;        // Number of points in the time grid.
    double _r;      // Interest rate
    double _sigma;  // Volatility
    double _T;      // Final time
    double _K;      // Strike Price
    double _R;      // Artificial limit of Asset price.

    void setSpaceGrid();
    void setTimeGrid();

    Vector* _s_grid;
    Vector* _t_grid;

    double _ds;
    double _dt;

};


#endif // ABSTRACTBLACKSCHOLES_HPP
