#include "PDE.h"
#include <cmath>




double PDE::getTheta() const
{
    double theta;
    switch (param.type)
    {
    case EXP:
        theta=0.;
        break;
    case IMP:
        theta=1.;
        break;
    case C_N:
        theta=0.5;
        break;
    default:
        throw "PDEType is not defined ";
        break;
    }
    return theta ;
}

/*
====================================================================
                     PDE COMPLETE 
====================================================================
*/
double PDE_C::get_cdt_bord_l(double t){
    if(get_payoff() == nullptr)
        throw "Payoff is null !";
    const auto * call = dynamic_cast<Call*>(get_payoff());
    if(call)
        return  0.0;
    return payoff->getStrike() * exp(-diffusion.r * (param.T - t));
}
double PDE_C::get_cdt_bord_L(double t ){
    if(get_payoff() == nullptr)
        throw "Payoff is null !";
    const auto * call = dynamic_cast<Call*>(get_payoff());
    if(call)
        return  param.L-payoff->getStrike() * exp(-diffusion.r * (param.T - t));
    return 0.0;

}


void PDE_C::computePDEcoeffs(int i,std::vector<double> & out) const {
    const double theta=getTheta();
    const double dt = dt_();
    const double sigmaSquared = diffusion.sigma*diffusion.sigma;
    const double iSquared= (double)i*(double)i;
    out.clear();
    out.resize(6);
    //coeff matrice multiplie par notre variable d interet
    out[0]= 1.+theta*dt*(diffusion.r+iSquared*sigmaSquared);//diag 
    out[1]= -0.5*theta*dt*(sigmaSquared*iSquared+diffusion.r*i) ;//upper
    out[2]= -0.5*theta*dt*(sigmaSquared*iSquared-diffusion.r*i);//lower 
    //coeff matrice du second membre
    out[3]= 1.-(1-theta)*dt*(diffusion.r+iSquared*sigmaSquared);//diag
    out[4]= 0.5*(1-theta)*dt*(sigmaSquared*iSquared+diffusion.r*i);//upper
    out[5]= 0.5*(1-theta)*dt*(sigmaSquared*iSquared-diffusion.r*i);//lower
 }
// m represente le passe (c l'etat qu'on veut calculer regressivement)
void PDE_C::computeBoardVect(int m, double & first,double & last){
    std::vector<double> coeff0,coeffL;
    computePDEcoeffs(1,coeff0);
    computePDEcoeffs(param.N-1,coeffL);
    const double dt =dt_();
    first= (coeff0[5]*get_cdt_bord_l(param.T-dt*(double)m)-coeff0[2]*get_cdt_bord_l(param.T-dt*(double)(m+1)));
    last = (coeffL[4]*get_cdt_bord_L(param.T-dt*(double)m)-coeffL[1]*get_cdt_bord_L(param.T-dt*(double)(m+1)));

}
void PDE_C::computeStartingVect(std::vector<double> & out){
    out.clear();
 
    out.resize(getNBSteps_S()-1);
    
    for (int i=0;i<getNBSteps_S()-1;i++){
        double ds = ds_();
        out[i]=get_payoff()->operator()(param.l+((double)i+1.)*ds);

    }
}
/*
====================================================================
                PDE REDUITE 
====================================================================
*/
double PDE_R::get_cdt_bord_l(double t){
    if(get_payoff() == nullptr)
        throw "Payoff is null !";
    const auto * call = dynamic_cast<Call*>(get_payoff());
    if(call)
        return  0.0;
    const double sigmaSquared= diffusion.sigma*diffusion.sigma;
    return payoff->getStrike() * exp(-diffusion.r * t) -exp((diffusion.r-sigmaSquared*0.5)*(param.T-t)+param.l);
}
double PDE_R::get_cdt_bord_L(double t ){
    if(get_payoff() == nullptr)
        throw "Payoff is null !";
    const double sigmaSquared= diffusion.sigma*diffusion.sigma;
    const auto * call = dynamic_cast<Call*>(get_payoff());
    if(call)
        return  exp((diffusion.r-sigmaSquared*0.5)*(param.T-t)+param.L)-payoff->getStrike() * exp(-diffusion.r * t);
    return 0.0;

}





void PDE_R::computePDEcoeffs(int /*i*/,std::vector<double> & out) const{

    const double theta=getTheta();
    const double dt = dt_();
    const double ds = ds_();
    const double dsSquared =ds*ds;
    const double sigmaSquared = diffusion.sigma*diffusion.sigma;
    const double alpha = 0.5*sigmaSquared*dt/dsSquared;
    out.clear();
    out.resize(6);
    //coeff matrice multiplie par notre variable d interet
    out[0]= 1+2*alpha*theta;//diag 
    out[1]= -alpha*theta;//upper
    out[2]= -alpha*theta;//lower 
    //coeff matrice du second membre
    out[3]= 1-2*alpha*(1-theta);//diag
    out[4]= alpha*(1-theta);//upper
    out[5]= alpha*(1-theta);//lower

}
// m represente le futur (c l'etat qu'on veut calculer progressivement )
void PDE_R::computeBoardVect(int m, double & first,double & last){
    std::vector<double> coeff0,coeffL;
    computePDEcoeffs(1,coeff0);
    computePDEcoeffs(param.N-1,coeffL);
    const double dt =dt_();
    first= (coeff0[5]*get_cdt_bord_l(dt*(double)(m))-coeff0[2]*get_cdt_bord_l(dt*(double)(m+1)));
    last = (coeffL[4]*get_cdt_bord_L(dt*(double)(m))-coeffL[1]*get_cdt_bord_L(dt*(double)(m+1)));

}

void PDE_R::computeStartingVect(std::vector<double> & out){
    out.clear();
    out.resize(getNBSteps_S()-1);
    for (int i=0;i<getNBSteps_S()-1;i++){
        
        const double ds = ds_();
        const double sigmaSquared = diffusion.sigma*diffusion.sigma;
        const double s = exp((diffusion.r-sigmaSquared*0.5)*param.T+param.l +(double)(i+1)*ds);
        out[i] = payoff->operator()(s);

    }
}

void PDE_R::normalizeSolution(std::vector<double> & sol) const{
    const double normalizingCoeff = exp(-diffusion.r*param.T);
    for (double &element : sol) {
        element *= normalizingCoeff;
    }
}