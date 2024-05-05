#ifndef DISCRETISATION_H
#define DISCRETISATION_H
#include <algorithm>

#include <iostream>
#include <vector>
#include <cmath>
enum STYPE {STANDARD, LOG};
class Discretisation
{
private:
    double l,L;
    int N;
    std::vector<double> grid;
    STYPE type{STYPE::STANDARD};

public:
    Discretisation(double l, double L,int N,STYPE t);
    void configureGrid();
    std::vector<double> getGrid(){return grid;}
    void computeSolutionOnNewGrid(const std::vector<double> & sol,
                                  const std::vector<double> & newGrid,
                                  std::vector<double> & newSol);
};







#endif //Discretisation.h