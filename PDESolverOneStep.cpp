#include "PDESolverOneStep.h"

PDESolverOneStep::PDESolverOneStep(std::unique_ptr<PDE> pde):Mx(1),Mb(1),pde(std::move(pde))
{
    configureOneStep();
}
void PDESolverOneStep::configureOneStep(){
    if(pde==nullptr)
        throw "pde is NULL!";
    int n =pde->getNBSteps_S()-1;
    Mx=Tridiag(n);
    Mb=Tridiag(n);
    
    const auto * pde_r =dynamic_cast<PDE_R*>(pde.get());
    if (pde_r){
        std::vector<double> coeffs;
        int dummy(0);
        pde->computePDEcoeffs(dummy,coeffs);
        for(int i=0;i<n;++i){
            Mx.set(i,coeffs[0],coeffs[1],coeffs[2]);
            Mb.set(i,coeffs[3],coeffs[4],coeffs[5]);
        }
    }
    else{
        
        for(int i=0;i<n;++i){
            std::vector<double> coeffs_i;
            pde->computePDEcoeffs(i,coeffs_i);
            Mx.set(i,coeffs_i[0],coeffs_i[1],coeffs_i[2]);
            Mb.set(i,coeffs_i[3],coeffs_i[4],coeffs_i[5]);
        }
    }
}

std::vector<double> PDESolverOneStep::solveOneStep(const std::vector<double> & xprev ){
    if(xprev.size() != static_cast<std::size_t>(pde->getNBSteps_S() - 1))
        throw "the previous vector size is not N-1";
    // on calcule d'abord le second membre     
    std::vector<double> sec_memb=Mb.multiply(xprev);
    double first(0.),last(0.);
    pde->computeBoardVect(step,first,last);
    sec_memb[0]+=first;
    sec_memb[xprev.size()-1] += last;
    // puis la solution sera stocke dans sec_memb pour economiser la creation d'un nouvel std::vector
    sec_memb=Mx.solve(sec_memb);
    return sec_memb;

}
std::vector<double> PDESolver::solve(){
    std::vector<double> sol ;
    pde->computeStartingVect(sol);
    for (int i=0;i<pde->getNBSteps_T();++i){
        step=i;
        sol=solveOneStep(sol);

    }
    pde->normalizeSolution(sol);
    return sol;

}