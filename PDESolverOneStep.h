// PDESolverOneStep.h
#ifndef PDESOLVERONESTEP_H
#define PDESOLVERONESTEP_H

#include <cmath>
#include "PDE.h"
#include "Tridiag.h"


class PDESolverOneStep
{
protected: 
    Tridiag Mx;
    Tridiag Mb;
    std::unique_ptr<PDE> pde;
    mutable int step;//temps i+dt
    
    
public:
    PDESolverOneStep(std::unique_ptr<PDE> pde);
    void configureOneStep();
    int getstep() const {return step;}
    void setStep(int m)const {step=m;}
    std::vector<double> solveOneStep(const std::vector<double> & xprev ); 
    
    
};

class PDESolver: public PDESolverOneStep{
    public:
    using PDESolverOneStep::PDESolverOneStep;
    std::vector<double> solve();


};

#endif // PDESolverOneStep.h
