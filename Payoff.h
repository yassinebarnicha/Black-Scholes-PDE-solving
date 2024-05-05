// Payoff.h
#ifndef PAYOFF_H
#define PAYOFF_H

#include <algorithm> // Pour std::max

class Payoff {
protected:
 double k_ ;

public:
    
    Payoff(double k): k_(k){}
    virtual ~Payoff() {}
    virtual double operator()(double S) const = 0;
    double getStrike(){return k_;}
};

class Put : public Payoff {
public:
    Put(double k) : Payoff(k){}
    double operator()(double s) const override {
        return std::max(k_ - s, 0.);
    }

    
};

class Call : public Payoff {

public:
    Call(double k) : Payoff(k) {}
    double operator()(double s) const override {
        return std::max(s - k_, 0.);
    }

};

#endif //PAYOFF_H
