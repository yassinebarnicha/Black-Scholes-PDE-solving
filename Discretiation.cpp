#include "Discretisation.h"


Discretisation::Discretisation(double l, double L,int N,STYPE t):l(l),L(L),N(N),type(t)
{
    configureGrid();
}

void Discretisation::configureGrid(){
    grid.resize(N-1);
    double ds(0.0);
    switch (type)
    {
    case STANDARD:
        ds = (L-l)/N;
        for (int i=0;i<N-1;++i){
            grid[i]=l+ds*(double)(i+1);
        }
        break;
    case LOG:
        ds = log(L/l)/N;
        for (int i=0;i<N-1;++i){
            grid[i]=exp(log(l)+ds*(double)(i+1));
        }



        break;
    default:
        throw "type is not handled !";
    }
}


void Discretisation::computeSolutionOnNewGrid(const std::vector<double> & sol,const std::vector<double> & newGrid,std::vector<double> & newSol)
{   
    newSol.clear();
    newSol.resize(sol.size());
    for (size_t i = 0; i < newGrid.size(); ++i) {
        // Trouver l'index de la valeur la plus proche dans oldGrid
        auto lower = std::lower_bound(grid.begin(), grid.end(), newGrid[i]);
        size_t index = std::distance(grid.begin(), lower);

        if (index == 0) {
            // Bord inférieur de la grille
            newSol[i] = sol.front();
        } else if (index >= grid.size()) {
            // Bord supérieur de la grille
            newSol[i] = sol.back();
        } else {
            // Interpolation linéaire
            double t = (newGrid[i] - grid[index - 1]) / (grid[index] - grid[index - 1]);
            newSol[i] = sol[index - 1] + t * (sol[index] - sol[index - 1]);
        }
    }

}


