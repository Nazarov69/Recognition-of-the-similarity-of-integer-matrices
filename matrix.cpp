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
        error(-1);
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
        nod = NOD(characteristic_matrix_A[2], characteristic_matrix_A[3]);
        eigenvector[0] = characteristic_matrix_A[2];
        eigenvector[1] = characteristic_matrix_A[3];
    }
    else {
        nod = 1;
        eigenvector[0] = 0;
        eigenvector[1] = 1;
        /*S = { 1,0,0,1 };
        std::vector<int> inverse_S(dimension * dimension, 0);
        std::vector<int> left_operation(dimension * dimension, 0);
        std::vector<int> C(dimension * dimension, 0);
        inverse_S = inverse_matrix(S);
        left_operation = multiplication(inverse_S, matr);
        C = multiplication(left_operation, S);

        printf("input matrix:\n");
        for (int i = 0; i < dimension; i++) {
            printf("| %d %d |\n", matr[i * dimension], matr[i * dimension + 1]);
        }

        if (eigenvalue_A1 == eigenvalue_A2) {
            printf("\na = %d\n", eigenvalue);
        }
        else {
            printf("\na = %d, b = %d\n", eigenvalue_A1, eigenvalue_A2);
        }

        printf("\nS^(-1) * A * S = C\n");
        for (int i = 0; i < dimension; i++) {
            printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
                inverse_S[i * dimension], inverse_S[i * dimension + 1],
                matr[i * dimension], matr[i * dimension + 1],
                S[i * dimension], S[i * dimension + 1],
                C[i * dimension], C[i * dimension + 1]);
        }
        printf("\n----------------------------------------------------------------\n\n");
        return C;*/
    }
    while (nod != 1) {
        eigenvector[0] /= nod;
        eigenvector[1] /= nod;
        nod = NOD(eigenvector[0], eigenvector[1]);
    }

    S[0] = eigenvector[1];
    S[2] = -eigenvector[0];
    int a = S[0];
    int b = S[2];
    int x, y;
    gcd(a, b, x, y);
    S[1] = -y;
    S[3] = x;
    parts_transform_matrix.push_back(S);
    
    std::vector<int> inverse_S(dimension * dimension, 0);
    std::vector<int> left_operation(dimension * dimension, 0);
    std::vector<int> C(dimension * dimension, 0);
    inverse_S = inverse_matrix(S);
    left_operation = multiplication(inverse_S, matr);
    C = multiplication(left_operation, S);

    printf("input matrix:\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d |\n", matr[i * dimension], matr[i * dimension + 1]);
    }

    if (eigenvalue_A1 == eigenvalue_A2) {
        printf("\na = %d\n", eigenvalue);
    }
    else {
        printf("\na = %d, b = %d\n", eigenvalue_A1, eigenvalue_A2);
    }

    printf("\nS^(-1) * A * S = C\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            inverse_S[i * dimension], inverse_S[i * dimension + 1],
            matr[i * dimension], matr[i * dimension + 1],
            S[i * dimension], S[i * dimension + 1],
            C[i * dimension], C[i * dimension + 1]);
    }

    if (C[1] < 0) { similar_matrix_with_negative_element(C); }

    if (eigenvalue_A1 != eigenvalue_A2 && C[1] >= C[3] - C[0]) {
        first_interval(C);
    }

    if (eigenvalue_A1 != eigenvalue_A2 && C[1] > (C[3] - C[0]) / 2) {
        second_interval(C);
    }
    printf("\n----------------------------------------------------------------\n");
    return C;
}

void similar::similar_matrix_with_negative_element(std::vector<int>& matr) {
    std::vector<int> T(dimension * dimension);
    T[0] = -1;
    T[1] = 0;
    T[2] = 0;
    T[3] = 1;
    parts_transform_matrix.push_back(T);
    std::vector<int> first_op(dimension * dimension, 0);
    std::vector<int> second_op(dimension * dimension, 0);
    std::vector<int> inverse_T(dimension * dimension);
    inverse_T = inverse_matrix(T);
    first_op = multiplication(inverse_T, matr);
    second_op = multiplication(first_op, T);
    matr = second_op;
    printf("\nsimilar matrix A:\n");
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
    parts_transform_matrix.push_back(T);
    std::vector<int> inverse_T(dimension * dimension);
    std::vector<int> first_operation(dimension * dimension, 0);
    std::vector<int> second_operation(dimension * dimension, 0);
    inverse_T = inverse_matrix(T);
    first_operation = multiplication(inverse_T, matr);
    second_operation = multiplication(first_operation, T);

    printf("\nfirst interval: T^(-1) * C * T = \n");
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
    parts_transform_matrix.push_back(R);
    std::vector<int> inverse_R(dimension * dimension);
    std::vector<int> first_operation(dimension * dimension, 0);
    std::vector<int> second_operation(dimension * dimension, 0);
    inverse_R = inverse_matrix(R);
    first_operation = multiplication(inverse_R, matr);
    second_operation = multiplication(first_operation, R);

    printf("\nsecond interval: R^(-1) * C * R = \n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            R[i * dimension], R[i * dimension + 1],
            matr[i * dimension], matr[i * dimension + 1],
            R[i * dimension], R[i * dimension + 1],
            second_operation[i * dimension], second_operation[i * dimension + 1]);
    }

    matr = second_operation;
}

void similar::class_comparison(const std::vector<int>& a, const std::vector<int>& c, const std::vector<int>& b, const std::vector<int>& d, const similar& sim) {
    if (c != d) {
        printf("\ndifferent classes\n");
        printf("\n----------------------------------------------------------------\n");
        return;
    }
    printf("\nidentical classes\n");
    printf("\n----------------------------------------------------------------\n");

    //******************************************************************************************************************

    std::vector<std::vector<int>> transform_matrix_A(parts_transform_matrix.size());
    std::vector<std::vector<int>> transform_matrix_B(sim.parts_transform_matrix.size());

    transform_matrix_A = parts_transform_matrix;
    transform_matrix_B = sim.parts_transform_matrix;

    std::vector<int> tr_A = transform_matrix_A[0];
    std::vector<int> tr_B = transform_matrix_B[0];

    for (int i = 0; i < transform_matrix_A.size() - 1; i++) {
        tr_A = { 0,0,0,0 };
        tr_A = multiplication(transform_matrix_A[i], transform_matrix_A[i + 1]);
        transform_matrix_A[i + 1] = tr_A;
    }

    for (int i = 0; i < transform_matrix_B.size() - 1; i++) {
        tr_B = { 0,0,0,0 };
        tr_B = multiplication(transform_matrix_B[i], transform_matrix_B[i + 1]);
        transform_matrix_B[i + 1] = tr_B;
    }

    printf("\nfirst transform matrix:\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d |\n", tr_A[i * dimension], tr_A[i * dimension + 1]);
    }

    printf("\nsecond transform matrix:\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d |\n", tr_B[i * dimension], tr_B[i * dimension + 1]);
    }

    //******************************************************************************************************************

    std::vector<int> inv(dimension * dimension);
    for (int i = sim.parts_transform_matrix.size() - 1; i >= 0; i--) {
        inv = inverse_matrix(sim.parts_transform_matrix[i]);
        parts_transform_matrix.push_back(inv);
    }

    std::vector<int> transform_matrix(dimension * dimension, 0);
    for (int i = 0; i < parts_transform_matrix.size() - 1; i++) {
        transform_matrix = { 0, 0, 0, 0 };
        transform_matrix = multiplication(parts_transform_matrix[i], parts_transform_matrix[i + 1]);
        parts_transform_matrix[i + 1] = transform_matrix;
    }

    std::vector<int> inverse(dimension * dimension, 0);
    std::vector<int> first_op(dimension * dimension, 0);
    std::vector<int> second_op(dimension * dimension, 0);
    inverse = inverse_matrix(transform_matrix);
    first_op = multiplication(inverse, a);
    second_op = multiplication(first_op, transform_matrix);
    printf("\ntransform matrix:\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d |\n", transform_matrix[i * dimension], transform_matrix[i * dimension + 1]);
    }

    //******************************************************************************************************************
    std::vector<int> inv1(dimension * dimension, 0);
    std::vector<int> f1_op(dimension * dimension, 0);
    std::vector<int> sec1_op(dimension * dimension, 0);
    inv1 = inverse_matrix(tr_A);
    f1_op = multiplication(inv1, a);
    sec1_op = multiplication(f1_op, tr_A);

    std::vector<int> inv2(dimension * dimension, 0);
    std::vector<int> f2_op(dimension * dimension, 0);
    std::vector<int> sec2_op(dimension * dimension, 0);
    inv2 = inverse_matrix(tr_B);
    f2_op = multiplication(inv2, b);
    sec2_op = multiplication(f2_op, tr_B);
    
    printf("\nmatrix check:\n");

    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            inv1[i * dimension], inv1[i * dimension + 1],
            a[i * dimension], a[i * dimension + 1],
            tr_A[i * dimension], tr_A[i * dimension + 1],
            sec1_op[i * dimension], sec1_op[i * dimension + 1]);
    }

    printf("\n");
    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            inv2[i * dimension], inv2[i * dimension + 1],
            b[i * dimension], b[i * dimension + 1],
            tr_B[i * dimension], tr_B[i * dimension + 1],
            sec2_op[i * dimension], sec2_op[i * dimension + 1]);
    }
    printf("\n");
    //******************************************************************************************************************

    for (int i = 0; i < dimension; i++) {
        printf("| %d %d | \t | %d %d | \t | %d %d | \t = \t | %d %d |\n",
            inverse[i * dimension], inverse[i * dimension + 1],
            a[i * dimension], a[i * dimension + 1],
            transform_matrix[i * dimension], transform_matrix[i * dimension + 1],
            second_op[i * dimension], second_op[i * dimension + 1]);
    }

    printf("\n----------------------------------------------------------------\n");
    if (sec1_op == c) {
        printf("\n1.DONE\n");
    }
    else {
        printf("\n1.NOT DONE\n");
    }
    if (sec2_op == d) {
        printf("\n2.DONE\n");
    }
    else {
        printf("\n2.NOT DONE\n");
    }
    if (second_op == b) {
        printf("\n3.DONE\n");
    }
    else {
        printf("\n3.NOT DONE\n");
    }
}

std::vector<int> similar::inverse_matrix(const std::vector<int>& matr) {
    std::vector<int> inverse(dimension * dimension);
    int det = matr[0] * matr[3] - matr[1] * matr[2];
    if (det != 1 && det != -1) {
        error(-2);
    }
    inverse[0] = matr[3] * det;
    inverse[1] = -matr[1] * det;
    inverse[2] = -matr[2] * det;
    inverse[3] = matr[0] * det;
    return inverse;
}

std::vector<int> similar::multiplication(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> mult(dimension * dimension, 0);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            for (int k = 0; k < dimension; k++) {
                mult[i * dimension + j] += a[i * dimension + k] * b[k * dimension + j];
            }
        }
    }
    return mult;
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