// Tridiag.h
#ifndef TRIDIAG_H
#define TRIDIAG_H

#include <vector>
#include <stdexcept>
class Tridiag {
private:
    std::vector<double> diag, upper, lower; // Diagonale principale, supérieure et inférieure
    int n; // Taille de la matrice

public:
    // Constructeur
    Tridiag(int size) : diag(size), upper(size - 1), lower(size - 1), n(size) {}

    // Setter pour les éléments
    void set(int i, double diagValue, double upperValue, double lowerValue);

    // Résoudre Ax = b avec l'algorithme de Thomas
    std::vector<double> solve(const std::vector<double>& b) const;
    std::vector<double> multiply(const std::vector<double>& vec) const ;
};


#endif // Tridiag.h
