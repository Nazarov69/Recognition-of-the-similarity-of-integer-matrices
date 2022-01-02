#include <iostream>
#include <vector>

#include "matrix.h"

int main() {
    std::vector<int> A(dimension * dimension);
    std::vector<int> C(dimension * dimension);
    A[0] = 1;
    A[1] = 5;
    A[2] = 0;
    A[3] = 2;

    similar sim;
    C = sim.similar_matrix(A);
    
    if (C[1] < 0) { sim.similar_matrix_with_negative_element(C); }

    if (sim.get_eigenvalue_A1() != sim.get_eigenvalue_A2() && C[1] >= C[3] - C[0]) {
        sim.first_interval(C);
    }

    if (sim.get_eigenvalue_A1() != sim.get_eigenvalue_A2() && C[1] > (C[3] - C[0]) / 2) {
        sim.second_interval(C);
    }
}