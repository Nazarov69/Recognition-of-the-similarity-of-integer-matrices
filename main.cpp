#include "matrix.h"

int main() {
    std::vector<int> A(dimension * dimension);
    std::vector<int> C(dimension * dimension);
    A[0] = -5;
    A[1] = 12;
    A[2] = 0;
    A[3] = 2;

    similar sim;
    C = sim.similar_matrix(A);

    std::vector<int> B(dimension * dimension);
    std::vector<int> D(dimension * dimension);
    B[0] = -5;
    B[1] = 23;
    B[2] = 0;
    B[3] = 2;

    similar sim_2;
    D = sim_2.similar_matrix(B);

    sim.class_comparison(A, C, B, D, sim_2);
}