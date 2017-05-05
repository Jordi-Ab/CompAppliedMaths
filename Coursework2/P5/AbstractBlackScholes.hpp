#ifndef ABSTRACTBLACKSCHOLES_HPP
#define ABSTRACTBLACKSCHOLES_HPP

#include <cmath>
#include "Vector.hpp"
#include "Helpers.hpp"
#include <stdexcept>

class AbstractBlackScholes{
public:

    // GETTERS

    // Returns the number of space points in
    // the space grid.
    double getNumSpacePoints();

    // Returns the number of space points in
    // the time grid.
    double getNumTimePoints();

    // Returns the value of the interest rate r.
    double getInterestRate();

    // Returns the value of the volatility sigma.
    double getVolatility();

    // Returns the final time T
    double getFinalTime();

    // Returns the Strike Price K
    double getStrikePrice();

    /*
     * Returns the final space point R.
     * Space points are prices of Asset, which
     * theoretically range from 0 to infinity,
     * but practically an Asset Price Limit is used,
     * which is usually taken as 3 or 4 times the
     * Strike Price of the option.
     */
    double getAssetPriceLimit();

    // Returns a copy of the Grid Vector of space points.
    Vector getSpaceGrid();

    // Returns a copy of the Grid Vector of time points.
    Vector getTimeGrid();

    // Returns the value of the spacing between space points.
    // Space Grid is uniform
    double getDeltaS();

    // Returns the value of the spacing between time points.
    // Time Grid is uniform
    double getDeltaT();

    // SETTERS

    // Sets a new value for the number of points in the
    // space grid. Resets the space grid after doing this.
    void setNumSpacePoints(int Nx);

    // Sets a new value for the number of points in the
    // time grid. Resets the time grid after doing this.
    void setNumTimePoints(int Nt);

    // Sets a new value for the interest rate
    void setInterestRate(double r);

    // Sets a new value for the volatility
    void setVolatility(double sigma);

    // Sets a new value for the final Time T.
    // Resets the time grid after doing this.
    void setFinalTime(double T);

    // Sets a new value for the stike Price K.
    void setStrikePrice(double K);

    // Sets a new value for the last Space Point.
    // Resets the space grid after doing this.
    void setAssetPriceLimit(double R);

    // Saves the space grid on an output file
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

    /*
     * Modifies the space grid "_s_grid" instance,
     * to be a uniformly spaced grid between 0 and _R,
     * containing _Nx+1 points.
    */
    void setSpaceGrid();

     /*
      * Modifies the time grid "_t_grid" instance,
      * to be a uniformly spaced grid between 0 and _T,
      * containing _Nt+1 points.
     */
    void setTimeGrid();

    // Space grid instance.
    Vector* _s_grid;

    // Time grid instance.
    Vector* _t_grid;

    // Delta s. Spacing between each space point.
    double _ds;
    // Delta t. Spacing between each time point.
    double _dt;

};


#endif // ABSTRACTBLACKSCHOLES_HPP
