#include "Tridiag.h"
void Tridiag::set(int i, double diagValue, double upperValue, double lowerValue) 
{
    if (i >= n) 
        throw std::out_of_range("Index out of range");
    diag[i] = diagValue;
    if (i < n - 1)
         upper[i] = upperValue;
    if (i > 0)
         lower[i - 1] = lowerValue;
}




// Résoudre Ax = b avec l'algorithme de Thomas
std::vector<double> Tridiag::solve(const std::vector<double>& b) const 
{
    if (b.size() != static_cast<std::size_t>(n))
         throw std::invalid_argument("Invalid size for vector b");

    std::vector<double> c_prime(n), d_prime(n), x(n);
    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = b[0] / diag[0];

        // Forward sweep
    for (int i = 1; i < n; i++) {
        double m = 1.0 / (diag[i] - lower[i - 1] * c_prime[i - 1]);
        c_prime[i] = i < n - 1 ? upper[i] * m : 0;
        d_prime[i] = (b[i] - lower[i - 1] * d_prime[i - 1]) * m;
        
    }

    // Backward substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}



std::vector<double> Tridiag::multiply(const std::vector<double>& vec) const 
{
    if (vec.size() != static_cast<std::size_t>(n)) 
        throw std::invalid_argument("Invalid size for vector");

    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        result[i] += diag[i] * vec[i];
        if (i < n - 1) result[i] += upper[i] * vec[i + 1];
        if (i > 0) result[i] += lower[i - 1] * vec[i - 1];
    }

    return result;
}
