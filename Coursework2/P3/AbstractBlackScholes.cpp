#include "AbstractBlackScholes.hpp"

double AbstractBlackScholes::getNumSpacePoints(){return _Nx;}

double AbstractBlackScholes::getNumTimePoints(){return _Nt;}

double AbstractBlackScholes::getInterestRate(){return _r;}

double AbstractBlackScholes::getVolatility(){return _sigma;}

double AbstractBlackScholes::getFinalTime(){return _T;}

double AbstractBlackScholes::getStrikePrice(){return _K;}

double AbstractBlackScholes::getAssetPriceLimit(){return _R;}

Vector AbstractBlackScholes::getSpaceGrid(){
    Vector copy(*_s_grid);
    return copy;
}

Vector AbstractBlackScholes::getTimeGrid(){
    Vector copy(*_t_grid);
    return copy;
}

double AbstractBlackScholes::getDeltaS(){
    return _ds;
}

double AbstractBlackScholes::getDeltaT(){
    return _dt;
}

void AbstractBlackScholes::setNumSpacePoints(int Nx){
    _Nx = Nx;
    // Recompute the space grid.
    setSpaceGrid();
}

void AbstractBlackScholes::setNumTimePoints(int Nt){
    _Nt = Nt;
    // Recompute the time grid.
    setTimeGrid();
}

void AbstractBlackScholes::setInterestRate(double r){
    _r = r;
}

void AbstractBlackScholes::setVolatility(double sigma){
    _sigma = sigma;
}

void AbstractBlackScholes::setFinalTime(double T){
    _T = T;
    // Recompute the time grid.
    setTimeGrid();
}

void AbstractBlackScholes::setStrikePrice(double K){
    _K = K;
}

void AbstractBlackScholes::setAssetPriceLimit(double R){
    _R = R;
    // Recompute the space grid.
    setSpaceGrid();
}

void AbstractBlackScholes::setSpaceGrid(){
    Vector a_grid = hlp.linspace(0, _R, _Nx+1);
    _s_grid = new Vector(a_grid);
    _ds = a_grid[1] - a_grid[0];
}

void AbstractBlackScholes::setTimeGrid(){
    Vector a_grid = hlp.linspace(0, _T, _Nt+1);
    _t_grid = new Vector(a_grid);
    _dt = a_grid[1] - a_grid[0];
}

void AbstractBlackScholes::saveGrid(std::string file_name){
    std::ofstream out_file;
    hlp.openOutputFile(file_name, out_file);
    hlp.saveVectorToFile(*_s_grid, file_name);
    hlp.closeOutputFile(out_file);
}
