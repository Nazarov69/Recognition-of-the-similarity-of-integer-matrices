#include <iostream>
#include <vector>

int NOD(int a, int b) {
    while (b != 0) {
        int c = a % b;
        a = b;
        b = c;
    }
    return abs(a);
}

int main() {
    int dimension = 2;
    std::vector<int> A(dimension * dimension);
    A[0] = 5;
    A[1] = 3;
    A[2] = -3;
    A[3] = -1;
    std::vector<int> S(dimension * dimension);
    std::vector<int> characteristic_matrix_A(dimension * dimension);
    std::vector<int> characteristic_matrix_A_2(dimension * dimension);
    std::vector<int> eigenvector(dimension * dimension);
    int eigenvalue_А1;
    int eigenvalue_А2;

    std::vector<int> result(dimension);
    int polynomial_A2 = -A[0] - A[3];
    int polynomial_A3 = A[0] * A[3] - A[1] * A[2];

    int determinant = sqrt((pow(polynomial_A2, 2) - 4 * polynomial_A3));
    if (pow(determinant, 2) == pow(polynomial_A2, 2) - 4 * polynomial_A3) {
        eigenvalue_А1 = (-polynomial_A2 + determinant) / 2;
        eigenvalue_А2 = (-polynomial_A2 - determinant) / 2;
    }
    else {
        return -1;
    }

    characteristic_matrix_A = A;
    for (int i = 0; i < dimension; i++) {
        characteristic_matrix_A[i * dimension + i] -= eigenvalue_А1;
    }

    int nod;
    if (characteristic_matrix_A[0] != 0 || characteristic_matrix_A[1] != 0) {
        nod = NOD(characteristic_matrix_A[0], characteristic_matrix_A[1]);
        result[0] = characteristic_matrix_A[0];
        result[1] = characteristic_matrix_A[1];
    }
    else if (characteristic_matrix_A[2] != 0 || characteristic_matrix_A[3] != 0) {
        nod = NOD(characteristic_matrix_A[3], characteristic_matrix_A[3]);
        result[0] = characteristic_matrix_A[2];
        result[1] = characteristic_matrix_A[3];
    }
    while (nod != 1) {
        result[0] /= nod;
        result[1] /= nod;
        nod = NOD(result[0], result[1]);
    }

    eigenvector[0] = -result[1];
    eigenvector[2] = result[0];

    if (eigenvector[0] < 0 && eigenvector[2] < 0) {
        eigenvector[0] *= -1;
        eigenvector[2] *= -1;
    }

    int a = std::abs(eigenvector[0]);
    int b = std::abs(eigenvector[2]);
    int p = 1, q = 0, r = 0, s = 1, x, y;
    while (a && b) {
        if (a >= b) {
            a = a - b;
            p = p - r;
            q = q - s;
        }
        else {
            b = b - a;
            r = r - p;
            s = s - q;
        }
    }
    if (a) {
        eigenvector[1] = q;
        eigenvector[3] = std::abs(p);
    }
    else {
        eigenvector[1] = s;
        eigenvector[3] = std::abs(r);
    }

    S = eigenvector;

    int Det = S[0] * S[3] - S[1] * S[2];
    if (Det != 1 && Det != -1) {
        return -2;
    }

    std::vector<int> inverse_matrix(dimension * dimension);
    inverse_matrix[0] = S[3] * Det;
    inverse_matrix[1] = -S[1] * Det;
    inverse_matrix[2] = -S[2] * Det;
    inverse_matrix[3] = S[0] * Det;

    std::vector<int> similar_matrix(dimension * dimension, 0);
    std::vector<int> left(dimension * dimension, 0);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                left[i * dimension + j] += inverse_matrix[i * dimension + k] * A[k * dimension + j];
            }
        }
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                similar_matrix[i * dimension + j] += left[i * dimension + k] * S[k * dimension + j];
            }
        }
    }

    printf("matrix A:\n");
    for (int i = 0; i < dimension; i++) {
        printf("%d %d\n", A[i * dimension], A[i * dimension + 1]);
    }
    printf("\n");
    printf("similar matrix:\n");
    for (int i = 0; i < dimension; i++) {
        printf("%d %d\n", similar_matrix[i * dimension], similar_matrix[i * dimension + 1]);
    }

    if (eigenvalue_А1 != eigenvalue_А2) {
        std::vector<int> result_2(dimension);
        std::vector<int> eigenvector_2(dimension * dimension);
        characteristic_matrix_A_2 = A;
        for (int i = 0; i < dimension; i++) {
            characteristic_matrix_A_2[i * dimension + i] -= eigenvalue_А2;
        }
        int nod;
        if (characteristic_matrix_A_2[0] != 0 || characteristic_matrix_A_2[1] != 0) {
            nod = NOD(characteristic_matrix_A_2[0], characteristic_matrix_A_2[1]);
            result_2[0] = characteristic_matrix_A_2[0];
            result_2[1] = characteristic_matrix_A_2[1];
        }
        else if (characteristic_matrix_A_2[2] != 0 || characteristic_matrix_A_2[3] != 0) {
            nod = NOD(characteristic_matrix_A_2[3], characteristic_matrix_A_2[3]);
            result_2[0] = characteristic_matrix_A_2[2];
            result_2[1] = characteristic_matrix_A_2[3];
        }
        while (nod != 1) {
            result_2[0] /= nod;
            result_2[1] /= nod;
            nod = NOD(result_2[0], result_2[1]);
        }

        eigenvector_2[0] = -result_2[1];
        eigenvector_2[2] = result_2[0];

        if (eigenvector_2[0] < 0 && eigenvector_2[2] < 0) {
            eigenvector_2[0] *= -1;
            eigenvector_2[2] *= -1;
        }

        int a = std::abs(eigenvector_2[0]);
        int b = std::abs(eigenvector_2[2]);
        int p = 1, q = 0, r = 0, s = 1, x, y;
        while (a && b) {
            if (a >= b) {
                a = a - b;
                p = p - r;
                q = q - s;
            }
            else {
                b = b - a;
                r = r - p;
                s = s - q;
            }
        }
        if (a) {
            eigenvector_2[1] = q;
            eigenvector_2[3] = std::abs(p);
        }
        else {
            eigenvector_2[1] = s;
            eigenvector_2[3] = std::abs(r);
        }

        std::vector<int> S_2 = eigenvector_2;

        int Det_2 = S_2[0] * S_2[3] - S_2[1] * S_2[2];
        if (Det_2 != 1 && Det_2 != -1) {
            return -2;
        }

        std::vector<int> inverse_matrix_2(dimension * dimension);
        inverse_matrix_2[0] = S_2[3] * Det_2;
        inverse_matrix_2[1] = -S_2[1] * Det_2;
        inverse_matrix_2[2] = -S_2[2] * Det_2;
        inverse_matrix_2[3] = S_2[0] * Det_2;

        std::vector<int> similar_matrix_2(dimension * dimension, 0);
        std::vector<int> left_2(dimension * dimension, 0);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                for (int k = 0; k < dimension; k++) {
                    left_2[i * dimension + j] += inverse_matrix_2[i * dimension + k] * similar_matrix[k * dimension + j];
                }
            }
        }
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                for (int k = 0; k < dimension; k++) {
                    similar_matrix_2[i * dimension + j] += left_2[i * dimension + k] * S_2[k * dimension + j];
                }
            }
        }

        printf("\n");
        printf("similar matrix 2:\n");
        for (int i = 0; i < dimension; i++) {
            printf("%d %d\n", similar_matrix_2[i * dimension], similar_matrix_2[i * dimension + 1]);
        }
    }
}