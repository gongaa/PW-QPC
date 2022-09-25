#include <iostream>
#include "Encoder/Encoder_RM.hpp"
#include "Encoder/Encoder_RM_CSS.hpp"
using namespace std;

// The X stabilizer matrix S_X = RM(m, m-rx-1)
// The Z stabilizer matrix S_Z = RM(m, m-rz-1)
// In total {m\choose \leq m-rx-1} + {m\choose \leq m-rz-1}
// = 2*2^m - {m\choose \leq rx} - {m\choose \leq rz}
// stabilizer generators
// So this CSS-constructed Quantum Reed-Muller code
// encodes K = {m\choose \leq rx} - {m\choose \leq rz} - 2^m
// logical qubits
Encoder_RM_CSS::Encoder_RM_CSS(const int& m, const int& rx, const int& rz)
{
    if (rx + rz < m - 1) {
        cerr << "Error: must satisfy rx+rz>=m-1" << endl;
        return;
    }
    this->m = m;
    this->rx = rx;
    this->rz = rz;
    this->N = 1 << m;
    this->K = Encoder_RM::calculate_K(m, rx) + Encoder_RM::calculate_K(m, rz) - (1 << m); 
}

// the logical X's are the quotient
// RM(m, rz) / RM(m, m-rx-1)
bool Encoder_RM_CSS::is_logical_X(const int *X_N)
{
    return Encoder_RM::is_logical(X_N, m, rz, m-rx-1);
}

// the logical Z's are the quotient
// RM(m, rx) / RM(m, m-rz-1)
bool Encoder_RM_CSS::is_logical_Z(const int *X_N)
{
    return Encoder_RM::is_logical(X_N, m, rx, m-rz-1);
}