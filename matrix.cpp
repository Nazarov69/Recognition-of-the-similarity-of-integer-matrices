#include "matrix.h"
std::vector<int> similar::similar_matrix(std::vector<int>& matr) {
    std::vector<int> S(dimension * dimension);
    std::vector<int> characteristic_matrix_A(dimension * dimension);
    std::vector<int> eigenvector(dimension);

    int polynomial_A2 = -matr[0] - matr[3];
    int polynomial_A3 = matr[0] * matr[3] - matr[1] * matr[2];

    int determinant = sqrt((pow(polynomial_A2, 2) - 4 * polynomial_A3));
    if (pow(determinant, 2) == pow(polynomial_A2, 2) - 4 * polynomial_A3) {
        eigenvalue_A1 = (-polynomial_A2 - determinant) / 2;
        eigenvalue_A2 = (-polynomial_A2 + determinant) / 2;
    }
    else {
        throw -1;
    }
    eigenvalue = std::min(eigenvalue_A1, eigenvalue_A2);

    characteristic_matrix_A = matr;
    for (int i = 0; i < dimension; i++) {
        characteristic_matrix_A[i * dimension + i] -= eigenvalue;
    }

    int nod;
    if (characteristic_matrix_A[0] != 0 || characteristic_matrix_A[1] != 0) {
        nod = NOD(characteristic_matrix_A[0], characteristic_matrix_A[1]);
        eigenvector[0] = characteristic_matrix_A[0];
        eigenvector[1] = characteristic_matrix_A[1];
    }
    else if (characteristic_matrix_A[2] != 0 || characteristic_matrix_A[3] != 0) {
        nod = NOD(characteristic_matrix_A[3], characteristic_matrix_A[3]);
        eigenvector[0] = characteristic_matrix_A[2];
        eigenvector[1] = characteristic_matrix_A[3];
    }
    while (nod != 1) {
        eigenvector[0] /= nod;
        eigenvector[1] /= nod;
        nod = NOD(eigenvector[0], eigenvector[1]);
    }

    S[0] = -eigenvector[1];
    S[2] = eigenvector[0];

    int a = std::abs(S[0]);
    int b = std::abs(S[2]);
    int x, y;
    gcd(a, b, x, y);
    S[1] = -y;
    S[3] = x;

    int Det = S[0] * S[3] - S[1] * S[2];
    if (Det != 1 && Det != -1) {
        throw -2;
    }

    std::vector<int> inverse_S(dimension * dimension);
    inverse_S[0] = S[3] * Det;
    inverse_S[1] = -S[1] * Det;
    inverse_S[2] = -S[2] * Det;
    inverse_S[3] = S[0] * Det;

    std::vector<int> C(dimension * dimension, 0);
    std::vector<int> left_operation(dimension * dimension, 0);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                left_operation[i * dimension + j] += inverse_S[i * dimension + k] * matr[k * dimension + j];
            }
        }
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                C[i * dimension + j] += left_operation[i * dimension + k] * S[k * dimension + j];
            }
        }
    }

    printf("matrix A:\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d |\n", matr[i * dimension], matr[i * dimension + 1]);
    }

    printf("\n");
    if (eigenvalue_A1 == eigenvalue_A2) {
        printf("a = %d\n", eigenvalue);
    }
    else {
        printf("a = %d, b = %d\n", eigenvalue_A1, eigenvalue_A2);
    }

    printf("\n");
    printf("S^(-1) * A * S = C\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            inverse_S[i * dimension], inverse_S[i * dimension + 1],
            matr[i * dimension], matr[i * dimension + 1],
            S[i * dimension], S[i * dimension + 1],
            C[i * dimension], C[i * dimension + 1]);
    }
    return C;
}

void similar::similar_matrix_with_negative_element(std::vector<int>& matr) {
    std::vector<int> T(dimension * dimension);
    T[0] = -1;
    T[1] = 0;
    T[2] = 0;
    T[3] = 1;

    std::vector<int> first_op(dimension * dimension, 0);
    std::vector<int> second_op(dimension * dimension, 0);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                first_op[i * dimension + j] += T[i * dimension + k] * matr[k * dimension + j];
            }
        }
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                second_op[i * dimension + j] += first_op[i * dimension + k] * T[k * dimension + j];
            }
        }
    }
    matr = second_op;
    printf("\n");
    printf("similar matrix A:\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d |\n", matr[i * dimension], matr[i * dimension + 1]);
    }
}

void similar::first_interval(std::vector<int>& matr) {
    std::vector<int> T(dimension * dimension, 0);
    T[0] = 1;
    T[1] = matr[1] / (matr[3] - matr[0]);
    T[2] = 0;
    T[3] = 1;

    std::vector<int> inverse_T(dimension * dimension);
    int d = T[0] * T[3] - T[1] * T[2];
    inverse_T[0] = T[3] * d;
    inverse_T[1] = -T[1] * d;
    inverse_T[2] = -T[2] * d;
    inverse_T[3] = T[0] * d;

    std::vector<int> first_operation(dimension * dimension, 0);
    std::vector<int> second_operation(dimension * dimension, 0);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                first_operation[i * dimension + j] += inverse_T[i * dimension + k] * matr[k * dimension + j];
            }
        }
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                second_operation[i * dimension + j] += first_operation[i * dimension + k] * T[k * dimension + j];
            }
        }
    }

    printf("\n");
    printf("first interval: T^(-1) * C * T = \n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            inverse_T[i * dimension], inverse_T[i * dimension + 1],
            matr[i * dimension], matr[i * dimension + 1],
            T[i * dimension], T[i * dimension + 1],
            second_operation[i * dimension], second_operation[i * dimension + 1]);
    }

    matr = second_operation;
}

void similar::second_interval(std::vector<int>& matr) {
    std::vector<int> R(dimension * dimension);
    R[0] = 1;
    R[1] = -1;
    R[2] = 0;
    R[3] = -1;

    std::vector<int> first_operation(dimension * dimension, 0);
    std::vector<int> second_operation(dimension * dimension, 0);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                first_operation[i * dimension + j] += R[i * dimension + k] * matr[k * dimension + j];
            }
        }
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                second_operation[i * dimension + j] += first_operation[i * dimension + k] * R[k * dimension + j];
            }
        }
    }

    printf("\n");
    printf("second interval: R^(-1) * C * R = \n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            R[i * dimension], R[i * dimension + 1],
            matr[i * dimension], matr[i * dimension + 1],
            R[i * dimension], R[i * dimension + 1],
            second_operation[i * dimension], second_operation[i * dimension + 1]);
    }
    matr = second_operation;
}

int similar::NOD(int a, int b) {
    while (b != 0) {
        int c = a % b;
        a = b;
        b = c;
    }
    return abs(a);
}

int similar::gcd(int a, int b, int& x, int& y) {
    if (a == 0) {
        x = 0; y = 1;
        return b;
    }
    int x1, y1;
    int d = gcd(b % a, a, x1, y1);
    x = y1 - b / a * x1;
    y = x1;
    return d;
}