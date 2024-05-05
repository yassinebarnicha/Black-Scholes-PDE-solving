#ifndef PDE_H
#define PDE_H
#include "Payoff.h"
#include <memory>
#include <vector>
#include <iostream>
struct Diffusion{
    double sigma, r;
};
enum PDEType {IMP,EXP,C_N};
struct PDEparam{
    int N,M;
    double T,L,l;
    PDEType type{PDEType::C_N};
};

class PDE
{
protected:
    std::unique_ptr<Payoff> payoff;
    Diffusion diffusion;
    PDEparam param;

public:
    PDE( std::unique_ptr<Payoff> payoff,Diffusion d,PDEparam pdeparam):payoff(std::move(payoff)),diffusion(d),param(pdeparam){}
    Payoff * get_payoff(){
        return payoff.get();
    }
    int getNBSteps_S(){return param.N;}
    int getNBSteps_T(){return param.M;}
    virtual double get_cdt_bord_l(double t) =0;
    virtual double get_cdt_bord_L(double t ) =0;
    double dt_() const {return param.T/((double)param.M);}
    double ds_()const {return (param.L-param.l)/((double)param.N);}
    double getTheta() const;
    virtual void computePDEcoeffs(int i,std::vector<double> & out) const =0;
    virtual void computeBoardVect(int m, double & first,double & last) =0;
    virtual void computeStartingVect(std::vector<double> & out) =0;
    virtual void normalizeSolution(std::vector<double> & /* sol */) const  {return;}
};



class PDE_C: public PDE{

 public:
    PDE_C( std::unique_ptr<Payoff> payoff,Diffusion d,PDEparam pdeparam):PDE(std::move(payoff),d,pdeparam){}
    double get_cdt_bord_l(double t) override;
    double get_cdt_bord_L(double t ) override;
    void computePDEcoeffs(int i,std::vector<double> & out) const override;
    void computeBoardVect(int m, double & first,double & last) override;
    void computeStartingVect(std::vector<double> & out) override;
    
};

class PDE_R: public PDE{
 public:
    PDE_R( std::unique_ptr<Payoff> payoff,Diffusion d,PDEparam pdeparam):PDE(std::move(payoff),d,pdeparam){}
    double get_cdt_bord_l(double t) override;
    double get_cdt_bord_L(double t ) override;
    void computePDEcoeffs(int i,std::vector<double> & out) const override;
    void computeBoardVect(int m, double & first,double & last) override;
    void computeStartingVect(std::vector<double> & out) override;
    void normalizeSolution(std::vector<double> & sol) const override;
    
};

#endif // PDE_H