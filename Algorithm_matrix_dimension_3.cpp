#include "Algorithm_matrix_dimension_3.h"

std::vector<std::vector<int>> Algorithm::algorithm(const std::vector<std::vector<int>>& matr, std::vector<std::vector<int>>& analog) {
    first_coefficient = 1;
    second_coefficient = -(matr[0][0] + matr[1][1] + matr[2][2]);
    third_coefficient = (get_det_2_2(matr[0][0], matr[0][1], matr[1][0], matr[1][1]) + 
        get_det_2_2(matr[1][1], matr[1][2], matr[2][1], matr[2][2]) + 
        get_det_2_2(matr[0][0], matr[0][2], matr[2][0], matr[2][2]));
    fourth_coefficient = -get_det_3_3(matr);

    get_eigenvalues(matr);

    std::vector<std::vector<int>> uni(DIM, std::vector<int>(DIM));
    std::vector<std::vector<int>> similar(DIM, std::vector<int>(DIM));
    std::vector<std::vector<int>> two_uni(DIM, std::vector<int>(DIM));
    std::vector<std::vector<int>> result(DIM, std::vector<int>(DIM));

    result = matr;

    uni = get_unimod_matrix(matr);

    similar = multiplication_and_inverse_B(matr, uni, analog);

    two_uni = get_two_uni(similar);

    result = multiplication_and_inverse_B(similar, two_uni, analog);

    if (first_eigenvalues == second_eigenvalues && first_eigenvalues == third_eigenvalues && second_eigenvalues == third_eigenvalues) {
        if (result[0][1] == 0 && result[0][2] == 0 && result[1][2] == 0) {
            printf("I\n a_1 = 0\n a_2 = 0\n a_3 = 0\n\n\n");
            return result;
        }
        else if (result[0][1] == 0 && (result[0][2] != 0 || result[1][2] != 0)) {
            int d = NOD(result[0][2], result[1][2]);
            std::vector<std::vector<int>> res(DIM, std::vector<int>(DIM));
            std::vector<std::vector<int>> S(DIM, std::vector<int>(DIM));
            S[0][0] = result[1][2] / d;
            S[0][1] = -result[0][2] / d;
            S[2][2] = 1;
            int u, v;
            gcd(result[0][2], result[1][2], u, v);
            if (get_det_2_2(S[0][0], S[0][1], u, v) == -1) {
                u = -u;
                v = -v;
                if (get_det_2_2(S[0][0], S[0][1], u, v) == -1) {
                    error("error u * a_2 + v * a_3 = d\n");
                }
            }
            S[1][0] = u;
            S[1][1] = v;
            std::vector<std::vector<int>> inverse_S(DIM, std::vector<int>(DIM));
            inverse_S = inverse_matrix(S);
            printf("a_1 = 0 and a_2 = 0 and a_3 = NOD(a_2, a_3)\n");
            res = multiplication_and_inverse_B(result, inverse_S, analog);
            if (res[0][1] != 0 || res[0][2] != 0 || res[1][2] != d) {
                error("error a_1 = 0 and a_2 = 0 and a_3 = NOD(a_2, a_3)\n");
            }
            printf("II\n a_1 = 0\n a_2 = 0\n a_3 >= 1\n\n\n");
            return res;
        }
        else if ((result[0][1] != 0 || result[0][2] != 0) && result[1][2] == 0) {
            int d = NOD(result[0][1], result[0][2]);
            std::vector<std::vector<int>> S(DIM, std::vector<int>(DIM));
            S[0][1] = 1;
            S[1][0] = result[0][2] / d;
            S[2][0] = -result[0][1] / d;
            int u, v;
            gcd(result[0][1], result[0][2], u, v);
            if (get_det_2_2(S[1][0], u, S[2][0], v) == -1) {
                u = -u;
                v = -v;
                if (get_det_2_2(S[1][0], u, S[2][0], v) == -1) {
                    error("error a_1 * u + a_2 * v = d\n");
                }
            }
            S[1][2] = u;
            S[2][2] = v;
            std::vector<std::vector<int>> res(DIM, std::vector<int>(DIM));
            printf("a_1 = 0 and a_2 = 0 and a_3 = NOD(a_1, a_2)\n");
            res = multiplication_and_inverse_B(result, S, analog);
            if (res[0][1] != 0 || res[0][2] != 0 || res[1][2] != d) {
                error("error a_1 = 0 and a_2 = 0 and a_3 = NOD(a_1, a_2)\n");
            }
            printf("II\n a_1 = 0\n a_2 = 0\n a_3 >= 1\n\n\n");
            return res;
        }
        else if (result[0][1] != 0 || result[1][2] != 0) {
            std::vector<std::vector<int>> res(DIM, std::vector<int>(DIM));
            std::vector<std::vector<int>> E(DIM, std::vector<int>(DIM));
            for (int i = 0; i < DIM; i++) {
                E[i][i] = 1;
            }
            if (result[0][1] < 0 && result[1][2] > 0) {
                E[0][0] = -1;
                printf("a_1 > 0\n");
                res = multiplication_and_inverse_B(result, E, analog);
                if (res[0][1] < 0) {
                    error("error a_1 > 0\n");
                }
            }
            else if (result[0][1] > 0 && result[1][2] < 0) {
                E[2][2] = -1;
                printf("a_3 > 0\n");
                res = multiplication_and_inverse_B(result, E, analog);
                if (res[1][2] < 0) {
                    error("error a_3 > 0\n");
                }
            }
            else if (result[0][1] < 0 && result[1][2] < 0) {
                E[1][1] = -1;
                printf("a_1 > and a_3 > 0\n");
                res = multiplication_and_inverse_B(result, E, analog);
                if (res[0][1] < 0 || res[1][2] < 0) {
                    error("error a_1 > 0 and a_3 > 0\n");
                }
            }
            else {
                res = result;
            }
            std::vector<std::vector<int>> S(DIM, std::vector<int>(DIM));
            for (int i = 0; i < DIM; i++) {
                S[i][i] = 1;
            }
            int nod = NOD(res[0][1], res[1][2]);
            int integer = res[0][2] < 0 && res[0][2] % nod != 0 ? (res[0][2] / nod) - 1 : res[0][2] / nod;
            int rem = res[0][2] - integer * nod;
            int q, e;
            gcd(res[0][1], res[1][2], e, q);
            q = -q;
            int det = get_det_2_2(res[0][1] * integer, q, res[1][2] * integer, e);
            if (det != rem - res[0][2] || det != -nod * integer) {
                q = -q;
                e = -e;
                det = get_det_2_2(res[0][1] * integer, q, res[1][2] * integer, e);
                if (det != rem - res[0][2] || det != -nod * integer) {
                    error("error a_1 * e - a_3 * q = r - a_2 = -ds\n");
                }
            }
            S[0][1] = q * integer;
            S[1][2] = e * integer;
            printf("0 <= a_2 < NOD(a_1, a_3)\n");
            res = multiplication_and_inverse_B(res, S, analog);
            if (0 > rem || rem >= nod || res[0][2] != rem) {
                error("error 0 <= a_2 < NOD(a_1, a_3)\n");
            }
            printf("III\n a_1 >= 1\n 0 <= a_2 < NOD(a_1, a_3) = %d\n a_3 >= 1\n\n\n", nod);
            return res;
        }
        return result;
    }
    if (pair(first_eigenvalues, second_eigenvalues, third_eigenvalues)) {
        std::vector<std::vector<int>> res(DIM, std::vector<int>(DIM));
        res = result;
        if (res[1][2] < 0) {
            std::vector<std::vector<int>> E_3(DIM, std::vector<int>(DIM));
            E_3[0][0] = 1;
            E_3[1][1] = 1;
            E_3[2][2] = -1;
            printf("a_3 > 0\n");
            res = multiplication_and_inverse_B(res, E_3, analog);
            if (res[1][2] < 0) {
                error("error a_3 > 0\n");
            }
        }
        int a_1 = res[0][1];
        int a_2 = res[0][2];
        int a_3 = res[1][2];

        int first_integer = get_integer(res[0][1], second_eigenvalues - first_eigenvalues);
        int first_rem = res[0][1] - first_integer * (second_eigenvalues - first_eigenvalues);

        if (res[0][1] >= abs(second_eigenvalues - first_eigenvalues) || res[0][1] < 0) {
            if (first_rem < 0 || first_rem >= abs(second_eigenvalues - first_eigenvalues)) {
                error("error 0 <= rem < |beta - alpha|\n");
            }
            std::vector<std::vector<int>> S(DIM, std::vector<int>(DIM));
            for (int i = 0; i < DIM; i++) {
                S[i][i] = 1;
            }
            S[0][1] = first_integer;
            printf("0 <= a_1 < |beta - alpha|\n");
            res = multiplication_and_inverse_B(res, S, analog);
            if (res[0][1] != first_rem || res[0][2] != a_2 - a_3 * first_integer) {
                error("error 0 <= a_1 < |beta - alpha|\n");
            }
        }
        if (res[0][1] < 0 || res[0][1] > (abs(second_eigenvalues - first_eigenvalues) / 2)) {
            /*a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];*/
            if (second_eigenvalues - first_eigenvalues > 0) {
                std::vector<std::vector<int>> T_1(DIM, std::vector<int>(DIM));
                T_1[0][0] = 1;
                T_1[0][1] = -1;
                T_1[1][1] = -1;
                T_1[2][2] = -1;
                printf("0 <= a_1 < |beta - alpha| / 2\n");
                res = multiplication_and_inverse_B(res, T_1, analog);
                if (res[0][1] != second_eigenvalues - first_eigenvalues - first_rem || res[0][2] != a_3 - a_2 + a_3 * first_integer) {
                    error("error 0 <= a_1 < |beta - alpha| / 2\n");
                }
            }
            if (second_eigenvalues - first_eigenvalues < 0) {
                std::vector<std::vector<int>> T_2(DIM, std::vector<int>(DIM));
                T_2[0][0] = 1;
                T_2[0][1] = 1;
                T_2[1][1] = -1;
                T_2[2][2] = -1;
                printf("0 <= a_1 < |beta - alpha| / 2\n");
                res = multiplication_and_inverse_B(res, T_2, analog);
                if (res[0][1] != first_eigenvalues - second_eigenvalues - first_rem || res[0][2] != -a_3 - a_2 + a_3 * first_integer) {
                    error("error 0 <= a_1 < |beta - alpha| / 2\n");
                }
            }
        }

        int nod = NOD(abs(second_eigenvalues - first_eigenvalues), res[0][1]);

        int second_integer = get_integer(res[0][2], nod);
        int second_rem = res[0][2] - second_integer * nod;
        if (second_rem < 0 || second_rem >= nod) {
            error("error 0 <= rem < nod\n");
        }
        int x_1, y_1;
        gcd(second_eigenvalues - first_eigenvalues, res[0][1], x_1, y_1);
        y_1 = -y_1;
        int det_1 = get_det_2_2(second_eigenvalues - first_eigenvalues, y_1, res[0][1], x_1);
        if (det_1 != nod) {
            x_1 = -x_1;
            y_1 = -y_1;
            det_1 = get_det_2_2(second_eigenvalues - first_eigenvalues, y_1, res[0][1], x_1);
            if (det_1 != nod) {
                error("error x_1 * (beta - alpha) - a_1 * y_1 = d\n");
            }
        }
        int x = second_integer * x_1;
        int y = second_integer * y_1;
        int det = get_det_2_2(second_eigenvalues - first_eigenvalues, y, res[0][1], x);
        if (det != nod * second_integer) {
            x = -x;
            y = -y;
            det = get_det_2_2((second_eigenvalues - first_eigenvalues) * second_integer, y, res[0][1] * second_integer, x);
            if (det != nod * second_integer) {
                error("error x * (beta - alpha) - a_1 * y = q * d\n");
            }
        }
        if (res[0][2] - det != second_rem) {
            error("error a_2 + a_1 * y - x * (deta - alpha) = rem\n");
        }

        a_1 = res[0][1];
        a_2 = res[0][2];
        a_3 = res[1][2];

        std::vector<std::vector<int>> T_3(DIM, std::vector<int>(DIM));
        for (int i = 0; i < DIM; i++) {
            T_3[i][i] = 1;
        }
        T_3[0][2] = x;
        T_3[1][2] = y;
        printf("a_2 = rem\n");
        res = multiplication_and_inverse_B(res, T_3, analog);
        if (res[0][2] != second_rem || res[0][1] != a_1 || res[0][2] != a_2 + a_1 * y - x * (second_eigenvalues - first_eigenvalues)) {
            error("error a_2 = rem\n");
        }

        if (abs(second_eigenvalues - first_eigenvalues) % 2 == 1 && res[0][1] != 0 && res[1][2] >= 1) {
            printf("II.a\n 1<= a_1 <= |beta - alpha| / 2 = %d\n 0 <= a_2 < d = %d\n a_3 >= 1\n\n\n",
                abs(second_eigenvalues - first_eigenvalues) / 2, nod);
            return res;
        }
        else if (abs(second_eigenvalues - first_eigenvalues) % 2 == 0 
            && res[0][1] != 0 && res[0][1] <= abs(second_eigenvalues - first_eigenvalues) / 2 - 1 &&  res[1][2] >= 1) {
            printf("II.b.1\n 1 <= a_1 <= |beta - alpha| / 2 - 1 = %d\n 0 <= a_2 < d = %d\n a_3 >= 1\n\n\n",
                abs(second_eigenvalues - first_eigenvalues) / 2 - 1, nod);
            return res;
        }

        if (res[0][1] == 0 && res[0][2] > abs(second_eigenvalues - first_eigenvalues) / 2) {
            std::vector<std::vector<int>> T_4(DIM, std::vector<int>(DIM));
            T_4[0][0] = -1;
            T_4[1][1] = 1;
            T_4[2][2] = 1;
            T_4[0][2] = x_1;
            printf("0 <= a_2 <= |beta - alpha| / 2\n");
            res = multiplication_and_inverse_B(res, T_4, analog);
            if (res[0][2] != nod - second_rem || res[0][2] != x_1 * (second_eigenvalues - first_eigenvalues) - second_rem) {
                error("error a_2 = nod - rem\n");
            }
            if (res[0][2] < 0 || res[0][2] > nod / 2 || res[0][2] > abs(second_eigenvalues - first_eigenvalues) / 2) {
                error("error 0 <= rem <= |beta - alpha| / 2\n");
            }
        }

        if (res[0][1] == 0 && 0 <= res[0][2] && res[0][2] <= abs(second_eigenvalues - first_eigenvalues) / 2 && res[1][2] >= 1) {
            printf("I\n a_1 = 0\n 0 <= a_2 <= |beta - alpha| / 2 = %d\n a_3 >= 1\n\n\n",
                abs(second_eigenvalues - first_eigenvalues) / 2);
            return res;
        }

        if (res[1][2] == 0 && (res[0][1] != 0 || res[0][2] != 0)) {
            int z = NOD(res[0][1], res[0][2]);
            std::vector<std::vector<int>> S_1(2, std::vector<int>(2));
            gcd(res[0][1], res[0][2], S_1[0][0], S_1[1][0]);
            int det = get_det_2_2(res[0][1], -S_1[1][0], res[0][2], S_1[0][0]);
            if (det != z) {
                S_1[0][0] = -S_1[0][0];
                S_1[1][0] = -S_1[1][0];
                det = get_det_2_2(res[0][1], -S_1[1][0], res[0][2], S_1[0][0]);
                if (det != z) {
                    error("error (a_1 a_2) * S_1 = (z 0)");
                }
            }
            S_1[0][1] = -res[0][2] / z;
            S_1[1][1] = res[0][2] == 0 ? 1 : -(res[0][1] * S_1[0][1]) / res[0][2];
            if (get_det_2_2(res[0][1], -S_1[1][1], res[0][2], S_1[0][1]) != 0) {
                error("error (a_1 a_2) * S_1 = (z 0)");
            }

            std::vector<std::vector<int>> T(DIM, std::vector<int>(DIM));
            T[0][0] = 1;
            for (int i = 1; i < DIM; i++) {
                for (int j = 1; j < DIM; j++) {
                    T[i][j] = S_1[i - 1][j - 1];
                }
            }
            printf("a_1 = NOD(a_1, a_2)\n");
            res = multiplication_and_inverse_B(res, T, analog);
            if (res[0][1] != z) {
                error("error a_1 = NOD(a_1, a_2)\n");
            }

            int d = NOD((second_eigenvalues - first_eigenvalues), z);
            int t_1, t_2;
            gcd((second_eigenvalues - first_eigenvalues), z, t_1, t_2);
            z = -z;
            if (get_det_2_2((second_eigenvalues - first_eigenvalues), z, t_2, t_1) != d) {
                t_1 = -t_1;
                t_2 = -t_2;
                if (get_det_2_2((second_eigenvalues - first_eigenvalues), z, t_2, t_1) != d) {
                    error("error (beta - alpha) * t_1 + z * t_2 = nod\n");
                }
            }
            std::vector<std::vector<int>> H(DIM, std::vector<int>(DIM));
            H[0][0] = 1;
            H[0][1] = -t_1;
            H[0][2] = z / d;
            H[1][1] = t_2;
            H[1][2] = (second_eigenvalues - first_eigenvalues) / d;
            if (NOD(H[1][1], H[1][2]) != 1) {
                error("error NOD(H[1][1], H[1][2]) = 1\n");
            }
            gcd(H[1][1], H[1][2], H[2][2], H[2][1]);
            H[1][2] = -H[1][2];
            if (get_det_2_2(H[1][1], H[1][2], H[2][1], H[2][2]) != 1) {
                H[2][1] = -H[2][1];
                H[2][2] = -H[2][2];
                if (get_det_2_2(H[1][1], H[1][2], H[2][1], H[2][2]) != 1) {
                    error("error H[1][1] * H[2][2] - H[1][2] * H[2][1]) = 1\n");
                }
            }
            printf("a_1 = NOD(beta - alpha, a_1, a_2)\n");
            res = multiplication_and_inverse_B(res, H, analog);
            if (res[0][1] != d) {
                error("error a_1 = NOD(beta - alpha, a_1, a_2)\n");
            }
        }

        if (res[1][2] == 0 && 
            (res[0][2] == 0 && res[0][1] == NOD((second_eigenvalues - first_eigenvalues),NOD(res[0][1], res[0][2])) || res[0][1] == 0)) {
            printf(" a_1 = d = %d\n a_2 = 0\n a_3 = 0\n\n\n",
                res[0][1] == 0 ? 0: NOD(second_eigenvalues - first_eigenvalues,NOD(res[0][1], res[0][2])));
            return res;
        }

        if (abs(second_eigenvalues - first_eigenvalues) % 2 == 0 && res[0][1] == abs(second_eigenvalues - first_eigenvalues) / 2) {
            int third_integer = get_integer(res[1][2], res[0][1]);
            int third_rem = res[1][2] - third_integer * res[0][1];
            if (second_eigenvalues - first_eigenvalues > 0) {
                if (third_rem / 2 + 1 <= res[0][2] && res[0][2] <= (res[0][1] + third_rem) / 2) {
                    int e = third_integer;
                    std::vector<std::vector<int>> S_1(DIM, std::vector<int>(DIM));
                    S_1[0][0] = -1;
                    S_1[0][1] = 1;
                    S_1[1][1] = 1;
                    S_1[1][2] = e;
                    S_1[2][2] = 1;
                    a_1 = res[0][1];
                    a_2 = res[0][2];
                    a_3 = res[1][2];
                    printf("-(a_1 - rem) / 2 <= a_2 <= rem / 2\n");
                    res = multiplication_and_inverse_B(res, S_1, analog);
                    if (res[0][1] != second_eigenvalues - first_eigenvalues - a_1 
                        || res[0][2] != a_2 - e * a_1 - a_2 || res[0][1] != a_1 || res[1][2] != a_3) {
                        error("error -(a_1 - rem) / 2 <= a_2 <= rem / 2\n");
                    }
                    if (res[0][2] != third_rem - a_2) {
                        error("error -(a_1 - rem) / 2 <= a_2 <= rem / 2\n");
                    }
                    if (-(res[0][1] - third_rem) / 2 > res[0][2] || res[0][2] > third_rem / 2) {
                        error("error -(a_1 - rem) / 2 <= a_2 <= rem / 2\n");
                    }
                }
                else if ((res[0][1] + third_rem) / 2 + 1 <= res[0][2] && res[0][2] <= res[0][1] - 1) {
                    std::vector<std::vector<int>> S_2(DIM, std::vector<int>(DIM));
                    S_2[0][0] = 1;
                    S_2[1][1] = 1;
                    S_2[1][2] = -1;
                    S_2[2][2] = 1;
                    a_1 = res[0][1];
                    a_2 = res[0][2];
                    a_3 = res[1][2];
                    printf("-(a_1 - rem) / 2 <= a_2 <= -1 <= rem / 2\n");
                    res = multiplication_and_inverse_B(res, S_2, analog);
                    if (res[0][2] != a_2 - a_1 || res[0][1] != a_1 || res[1][2] != a_3) {
                        error("error -(a_1 - rem) / 2 <= a_2 <= -1 <= rem / 2\n");
                    }
                    if (-(res[0][1] - third_rem) / 2 > res[0][2] || res[0][2] > third_rem / 2 || res[0][2] > -1) {
                        error("error -(a_1 - rem) / 2 <= a_2 <= -1 <= rem / 2\n");
                    }
                }
                printf("II.b.2\n b - a > 0\n a_1 = |beta - alpha| / 2 = %d\n %d = -(a_1 - r) / 2 <= a_2 <= r / 2 = %d , r = %d\n a_3 >= 1\n\n\n",
                    abs(second_eigenvalues - first_eigenvalues) / 2,
                    -(res[0][1] - third_rem) / 2, third_rem / 2, third_rem);
                return res;
            }
            else if (second_eigenvalues - first_eigenvalues < 0) {
                if ((res[0][1] - third_rem) / 2 + 1 <= res[0][2] <= res[0][1] - third_rem / 2 - 1) {
                    int e = third_integer + 1;
                    std::vector<std::vector<int>> S_3(DIM, std::vector<int>(DIM));
                    S_3[0][0] = 1;
                    S_3[0][1] = 1;
                    S_3[1][1] = -1;
                    S_3[1][2] = e;
                    S_3[2][2] = -1;
                    a_1 = res[0][1];
                    a_2 = res[0][2];
                    a_3 = res[1][2];
                    printf("-rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    res = multiplication_and_inverse_B(res, S_3, analog);
                    if (res[0][1] != first_eigenvalues - second_eigenvalues - a_1 || res[0][2] != e * a_1 - a_3 - a_2) {
                        error("error -rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    }
                    if (res[0][1] != a_1 || res[0][2] != e * a_1 - a_3 - a_2) {
                        error("error -rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    }
                    if (res[0][2] != a_1 - third_rem - a_2) {
                        error("error -rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    }
                    if (-(third_rem / 2) > res[0][2] || res[0][2] > (res[0][1] - third_rem) / 2) {
                        error("error -rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    }
                }
                else if (res[0][1] - third_rem / 2 <= res[0][2] && res[0][2] <= res[0][1] - 1) {
                    std::vector<std::vector<int>> S_4(DIM, std::vector<int>(DIM));
                    S_4[0][0] = 1;
                    S_4[1][1] = 1;
                    S_4[1][2] = -1;
                    S_4[2][2] = 1;
                    a_1 = res[0][1];
                    a_2 = res[0][2];
                    a_3 = res[1][2];
                    printf("-rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    res = multiplication_and_inverse_B(res, S_4, analog);
                    if (res[0][2] != a_2 - a_1) {
                        error("error -rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    }
                    if (-(third_rem / 2) > res[0][2] || res[0][2] > (res[0][1] - third_rem) / 2 || res[0][2] > -1) {
                        error("error -rem / 2 <= a_2 <= (a_1 - rem) / 2\n");
                    }
                }
                printf("II.b.2\n b - a < 0\n a_1 = |b - a| / 2 = %d\n %d = -r / 2 <= a_2 <= (a_1 - r) / 2 = %d , r = %d\n a_3 >= 1\n\n\n",
                    abs(second_eigenvalues - first_eigenvalues) / 2, -third_rem / 2,
                    (res[0][1] - third_rem) / 2, third_rem);
                return res;
            }
        }
        return res;
    }
    if (no_pair(first_eigenvalues, second_eigenvalues, third_eigenvalues)) {
        std::vector<std::vector<int>> res(DIM, std::vector<int>(DIM));
        res = result;
        int a_1 = res[0][1];
        int a_2 = res[0][2];
        int a_3 = res[1][2];

        int first_integer = get_integer(res[0][1], second_eigenvalues - first_eigenvalues);
        int first_rem = res[0][1] - first_integer * (second_eigenvalues - first_eigenvalues);

        if (res[0][1] >= second_eigenvalues - first_eigenvalues || res[0][1] < 0) {
            std::vector<std::vector<int>> S_1(DIM, std::vector<int>(DIM));
            for (int i = 0; i < DIM; i++) {
                S_1[i][i] = 1;
            }
            S_1[0][1] = first_integer;
            printf("0 <= a_1 < beta - alpha\n");
            res = multiplication_and_inverse_B(res, S_1, analog);
            if (res[0][1] != first_rem || res[0][2] != a_2 - a_3 * first_integer || res[1][2] != a_3) {
                error("error 0 <= a_1 < beta - alpha\n");
            }
        }

        if (res[0][1] > ((second_eigenvalues - first_eigenvalues) / 2) || res[0][1] < 0) {
            a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];
            std::vector<std::vector<int>> T_1(DIM, std::vector<int>(DIM));
            T_1[0][0] = 1;
            T_1[0][1] = -1;
            T_1[1][1] = -1;
            T_1[2][2] = -1;
            printf("0 <= a_1 <= (beta - alpha) / 2\n");
            res = multiplication_and_inverse_B(res, T_1, analog);
            if (res[0][1] != second_eigenvalues - first_eigenvalues - a_1 || res[0][2] != a_3 - a_2 || res[1][2] != a_3) {
                error("error 0 <= a_1 <= (beta - alpha) / 2\n");
            }
        }

        int third_integer = get_integer(res[1][2], third_eigenvalues - second_eigenvalues);
        int third_rem = res[1][2] - third_integer * (third_eigenvalues - second_eigenvalues);

        if (res[1][2] >= third_eigenvalues - second_eigenvalues || res[1][2] < 0) {
            a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];
            std::vector<std::vector<int>> S_2(DIM, std::vector<int>(DIM));
            for (int i = 0; i < DIM; i++) {
                S_2[i][i] = 1;
            }
            S_2[1][2] = third_integer;
            printf("0 <= a_3 < gamma - beta\n");
            res = multiplication_and_inverse_B(res, S_2, analog);
            if (res[0][1] != a_1 || res[0][2] != a_2 + a_1 * third_integer || res[1][2] != third_rem) {
                error("error 0 <= a_3 < gamma - beta\n");
            }
        }

        if (res[1][2] > ((third_eigenvalues - second_eigenvalues) / 2) || res[1][2] < 0) {
            a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];
            std::vector<std::vector<int>> T_2(DIM, std::vector<int>(DIM));
            T_2[0][0] = 1;
            T_2[1][1] = 1;
            T_2[1][2] = -1;
            T_2[2][2] = -1;
            printf("0 <= a_3 <= (gamma - beta) / 2\n");
            res = multiplication_and_inverse_B(res, T_2, analog);
            if (res[0][1] != a_1 || res[0][2] != -a_1 - a_2 || res[1][2] != third_eigenvalues - second_eigenvalues - a_3) {
                error("error 0 <= a_3 <= (gamma - beta) / 2\n");
            }
        }

        int second_integer = get_integer(res[0][2], third_eigenvalues - first_eigenvalues);
        int second_rem = res[0][2] - second_integer * (third_eigenvalues - first_eigenvalues);

        if (res[0][2] >= third_eigenvalues - first_eigenvalues || res[0][2] < 0) {
            a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];
            std::vector<std::vector<int>> S_3(DIM, std::vector<int>(DIM));
            for (int i = 0; i < DIM; i++) {
                S_3[i][i] = 1;
            }
            S_3[0][2] = second_integer;
            printf("0 <= a_2 < gamma - alpha\n");
            res = multiplication_and_inverse_B(res, S_3, analog);
            if (res[0][1] != a_1 || res[0][2] != second_rem || res[1][2] != a_3) {
                error("error 0 <= a_2 < gamma - alpha\n");
            }
        }

        if (res[0][1] != 0 && (third_eigenvalues - second_eigenvalues) % 2 == 1)
            printf("II.a\n 1 <= a_1 <= (beta - alpha) / 2 = %d\n 0 <= a_2 < gamma - alpha = %d\n 0 <= a_3 <= (gamma - beta) / 2 = %d\n\n\n",
                (second_eigenvalues - first_eigenvalues) / 2,
                third_eigenvalues - first_eigenvalues,
                (third_eigenvalues - second_eigenvalues) / 2);
        else if (res[0][1] != 0 && 
            (third_eigenvalues - second_eigenvalues) % 2 == 0 && 0 <= res[1][2] && res[1][2] <= (third_eigenvalues - second_eigenvalues) / 2 - 1)
            printf("II.b.1\n 1 <= a_1 <= (beta - alpha) / 2 = %d\n 0 <= a_2 < gamma - alpha = %d\n 0 <= a_3 <= (gamma - beta) / 2 - 1 = %d\n\n\n",
                (second_eigenvalues - first_eigenvalues) / 2,
                third_eigenvalues - first_eigenvalues,
                (third_eigenvalues - second_eigenvalues) / 2 - 1);

        if (res[0][1] == 0 && (res[0][2] > (third_eigenvalues - first_eigenvalues) / 2 || res[0][2] < 0)) {
            a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];
            std::vector<std::vector<int>> T_3(DIM, std::vector<int>(DIM));
            T_3[0][0] = 1;
            T_3[0][2] = -1;
            T_3[1][1] = -1;
            T_3[2][2] = -1;
            printf("a_1 = 0 and 0 <= a_2 <= (gamma - alpha) / 2\n");
            res = multiplication_and_inverse_B(res, T_3, analog);
            if (res[0][1] != 0 || res[0][2] != third_eigenvalues - first_eigenvalues - a_2 || res[1][2] != a_3) {
                error("error a_1 = 0 and 0 <= a_2 <= (gamma - alpha) / 2\n");
            }
        }

        if (res[0][1] == 0)
            printf("I\n a_1 = 0\n 0 <= a_2 <= (gamma - alpha) / 2 = %d\n 0 <= a_3 <= (gamma - beta) / 2 = %d\n\n\n",
                (third_eigenvalues - first_eigenvalues) / 2,
                (third_eigenvalues - second_eigenvalues) / 2);

        if ((third_eigenvalues - second_eigenvalues) % 2 == 0 && res[1][2] == (third_eigenvalues - second_eigenvalues) / 2 && res[0][1] != 0) {
            a_1 = res[0][1];
            a_2 = res[0][2];
            a_3 = res[1][2];
            if ((third_eigenvalues - first_eigenvalues - a_1) / 2 + 1 <= a_2 && a_2 <= third_eigenvalues - first_eigenvalues - a_1 / 2 - 1) {
                std::vector<std::vector<int>> T_4(DIM, std::vector<int>(DIM));
                T_4[0][0] = 1;
                T_4[0][2] = -1;
                T_4[1][1] = 1;
                T_4[1][2] = -1;
                T_4[2][2] = -1;
                printf("-a_1 / 2 <= a_2 <= (gamma - alpha - a_1) / 2\n");
                res = multiplication_and_inverse_B(res, T_4, analog);
                if (res[0][1] != a_1 || res[0][2] != third_eigenvalues - first_eigenvalues - a_1 - a_2 || 
                    res[1][2] != third_eigenvalues - second_eigenvalues - a_3 || res[1][2] != a_3) {
                    error("error - a_1 / 2 <= a_2 <= (gamma - alpha - a_1) / 2\n");
                }
            }
            else if (third_eigenvalues - first_eigenvalues - a_1 / 2 <= a_2 && a_2 <= third_eigenvalues - first_eigenvalues - 1) {
                std::vector<std::vector<int>> T_5(DIM, std::vector<int>(DIM));
                T_5[0][0] = 1;
                T_5[0][2] = 1;
                T_5[1][1] = 1;
                T_5[2][2] = 1;
                printf("-a_1 / 2 <= a_2 <= (gamma - alpha - a_1) / 2\n");
                res = multiplication_and_inverse_B(res, T_5, analog);
                if (res[0][1] != a_1 || res[0][2] != first_eigenvalues - third_eigenvalues + a_2 || res[1][2] != a_3) {
                    error("error - a_1 / 2 <= a_2 <= (gamma - alpha - a_1) / 2\n");
                }
            }
            printf("II.b.2\n 1 <= a_1 <= (beta - alpha) / 2 = %d\n %d = -a_1 / 2 <= a_2 <= (gamma - alpha - a_1) / 2 = %d\n a_3 = (gamma - beta) / 2 = %d\n\n\n",
                (second_eigenvalues - first_eigenvalues) / 2, -res[0][1] / 2,
                (third_eigenvalues - first_eigenvalues - res[0][1]) / 2,
                (third_eigenvalues - second_eigenvalues) / 2);
        }
        return res;
    }
    return result;
}

bool Algorithm::pair(int a, int b, int c) {
    if (a == b && a != c) return true;
    if (a == c && a != b) return true;
    if (b == c && b != a) return true;
    return false;
}

bool Algorithm::no_pair(int a, int b, int c) {
    if (a != b && a != c && b != c) return true;
    return false;
}

int Algorithm::get_det_2_2(const int& a, const int& b, const int& c,
    const int& d) {
    return a * d - b * c;
}

int Algorithm::get_det_3_3(const std::vector<std::vector<int>>& matr) {
    return matr[0][0] * get_det_2_2(matr[1][1], matr[1][2], matr[2][1], matr[2][2]) -
        matr[0][1] * get_det_2_2(matr[1][0], matr[1][2], matr[2][0], matr[2][2]) +
        matr[0][2] * get_det_2_2(matr[1][0], matr[1][1], matr[2][0], matr[2][1]);
}

void Algorithm::get_eigenvalues(const std::vector<std::vector<int>>& matr) {
    int** num = new int* [3];
    num[0] = &first_eigenvalues;
    num[1] = &second_eigenvalues;
    num[2] = &third_eigenvalues;

    int d = abs(fourth_coefficient);
    int i = 1;
    int k = 0;
    if (d == 0) {
        d = abs(third_coefficient);
        if (d == 0) {
            d = abs(second_coefficient);
            if (d == 0) {
                // error("input matrix 0\n");
                *num[0] = 0;
                *num[1] = 0;
                *num[2] = 0;
                k = 3;
            }
        }
    }
    while (i <= d && k < 3) {
        if (d % i == 0) {
            if (((first_coefficient * i + second_coefficient) * i + third_coefficient) * i + fourth_coefficient == 0) {
                *num[k++] = i;
            }
            if (((first_coefficient * (-i) + second_coefficient) * (-i) +third_coefficient) * (-i) + fourth_coefficient == 0) {
                *num[k++] = -i;
            }
        }
        i++;
    }
    if (k == 3) {
        if (*num[0] > *num[1]) {
            std::swap(*num[0], *num[1]);
        }
        if (*num[0] > *num[2]) {
            std::swap(*num[0], *num[2]);
        }
        if (*num[1] > *num[2]) {
            std::swap(*num[1], *num[2]);
        }
    }
    if (k == 1) {
        if (get_multiplicity(*num[0]) == 3) {
            *num[k++] = *num[0];
            *num[k++] = *num[0];
        }
        else {
            if (get_multiplicity(0) == 2) {
                *num[k++] = 0;
                *num[k++] = 0;
            }
            else {
                error("integers not found\n");
            }
        }
    }
    if (k == 2) {
        if (get_multiplicity(*num[0]) == 2) {
            *num[k++] = *num[0];
            std::swap(*num[0], *num[1]);
            if (*num[1] != *num[2]) {
                std::swap(*num[0], *num[1]);
            }
        }
        else if (get_multiplicity(*num[1]) == 2) {
            *num[k++] = *num[1];
            if (*num[1] != *num[2]) {
                std::swap(*num[0], *num[1]);
            }
        }
        else {
            if (get_multiplicity(0) == 1) {
                *num[k++] = 0;
                if (*num[0] > *num[1]) {
                    std::swap(*num[0], *num[1]);
                }
                if (*num[0] > *num[2]) {
                    std::swap(*num[0], *num[2]);
                }
                if (*num[1] > *num[2]) {
                    std::swap(*num[1], *num[2]);
                }
            }
            else {
                error("integers not found\n");
            }
        }
    }
    if (*num[0] == *num[1]) {
        std::swap(*num[1], *num[2]);
    }
    printf("lambda_1 = %d\nlambda_2 = %d\nlambda_3 = %d\n\n", first_eigenvalues, second_eigenvalues, third_eigenvalues);

    if (num != NULL) {
        delete[] num;
    }
}

int Algorithm::get_multiplicity(const int elem) {
    std::vector<int> coeffs;
    coeffs.push_back(first_coefficient);
    coeffs.push_back(second_coefficient);
    coeffs.push_back(third_coefficient);
    coeffs.push_back(fourth_coefficient);
    int multiplicity = 0;
    while (true) {
        for (int i = 1; i < coeffs.size(); i++) {
            coeffs[i] = elem * coeffs[i - 1] + coeffs[i];
        }
        if (coeffs[coeffs.size() - 1] == 0) {
            coeffs.erase(coeffs.begin() + coeffs.size() - 1);
            multiplicity++;
        }
        else {
            break;
        }
    }
    return multiplicity;
}

void Algorithm::print_coeffs() {
    printf("%d\n%d\n%d\n%d\n", first_coefficient, second_coefficient, third_coefficient, fourth_coefficient);
}

int Algorithm::NOD(int a, int b) {
    while (b != 0) {
        int c = a % b;
        a = b;
        b = c;
    }
    return abs(a);
}

int Algorithm::gcd(int a, int b, int& x, int& y) {
    int nod = NOD(a, b);
    while (nod != 1) {
        a /= nod;
        b /= nod;
        nod = NOD(a, b);
    }
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    int x1, y1;
    int d = gcd(b % a, a, x1, y1);
    x = y1 - b / a * x1;
    y = x1;
    return d;
}

std::vector<std::vector<int>> Algorithm::multiplication(const std::vector<std::vector<int>>& a, const std::vector<std::vector<int>>& b) {
    std::vector<std::vector<int>> mult(DIM, std::vector<int>(DIM));
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            for (int k = 0; k < DIM; k++) {
                mult[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return mult;
}

std::vector<std::vector<int>> Algorithm::multiplication_and_inverse_B(const std::vector<std::vector<int>>& a, 
    const std::vector<std::vector<int>>& b, std::vector<std::vector<int>>& analog) {
    std::vector<std::vector<int>> left_op;
    std::vector<std::vector<int>> right_op;
    std::vector<std::vector<int>> inverse_b = inverse_matrix(b);
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++)
            inverse_b[i][j] < 0 ? printf("%d ", inverse_b[i][j]) : printf("%d  ", inverse_b[i][j]);
        printf("|\t\t\t");
        for (int j = 0; j < DIM; j++)
            a[i][j] < 0 ? printf("%d ", a[i][j]) : printf("%d  ", a[i][j]);
        printf("|\t\t\t");
        for (int j = 0; j < DIM; j++)
            b[i][j] < 0 ? printf("%d ", b[i][j]) : printf("%d  ", b[i][j]);
        printf("\t|\t=\n");
    }
    left_op = multiplication(inverse_b, a);

    right_op = multiplication(left_op, b);
    analog = multiplication(analog, b);
    printf("similar matrix :\n");
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++)
            right_op[i][j] < 0 ? printf("%d ", right_op[i][j]) : printf("%d  ", right_op[i][j]);
        printf("\t|\t\t\n");
    }
    printf("\n");
    return right_op;
}

std::vector<std::vector<int>> Algorithm::inverse_matrix(const std::vector<std::vector<int>>& matr) {
    std::vector<std::vector<int>> inverse(DIM, std::vector<int>(DIM));
    inverse = matr;
    int det = get_det_3_3(matr);
    if (det != 1 && det != -1) {
        exit(777);
    }

    for (int i = 0, row, col; i < DIM; i++, row *= 2) {
        row = 1;
        col = 1;
        for (int j = 0; j < DIM; j++, col *= 2) {
            inverse[i][j] =
                pow(-1, i + j) *
                get_det_2_2(matr[i == 1 ? 0 : ((i + 1) % DIM)][j == 1 ? 0 : ((j + 1) % DIM)],
                    matr[i == 1 ? 0 : ((i + 1) % DIM)][j == 1 ? 2 : ((j + 2) % DIM)],
                    matr[i == 1 ? 2 : ((i + 2) % DIM)][j == 1 ? 0 : ((j + 1) % DIM)],
                    matr[i == 1 ? 2 : ((i + 2) % DIM)][j == 1 ? 2 : ((j + 2) % DIM)]);
        }
    }

    std::swap(inverse[1][0], inverse[0][1]);
    std::swap(inverse[2][0], inverse[0][2]);
    std::swap(inverse[2][1], inverse[1][2]);

    if (det == -1) {
        for (int i = 0, k; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                inverse[i][j] /= det;
            }
        }
    }
    return inverse;
}

std::vector<std::vector<int>> Algorithm::get_unimod_matrix(const std::vector<std::vector<int>>& matr) {
    std::vector<std::vector<int>> matrix_diff;
    std::vector<std::vector<int>> uni(DIM, std::vector<int>(DIM));
    matrix_diff = matr;

    int eigenvalue = first_eigenvalues;

    for (int i = 0; i < DIM; i++) {
        matrix_diff[i][i] -= eigenvalue;
    }
    std::vector<int> xx;

    xx = gauss(matrix_diff);

    int u, v;
    std::vector<int> two_col(DIM);

    int min = 0;
    if (abs(xx[min]) > abs(xx[1])) { 
        min = 1;
    }
    if (abs(xx[min]) > abs(xx[2])) {
        min = 2;
    }

    if (xx[1] == 0 && abs(xx[2]) == 1) {
        for (int i = 0; i < DIM; i++) {
            uni[i][0] = xx[i];
        }
        uni[0][1] = xx[2];
        uni[1][1] = 0;
        uni[2][1] = 0;
        uni[2][2] = 0;
        uni[0][2] = 0;
        uni[1][2] = 1;
    }
    else {
        while (true) {
            if (NOD(xx[(min + 1) % DIM], xx[(min + 2) % DIM]) != 1) {
                min = (min + 1) % DIM;
                if (NOD(xx[(min + 1) % DIM], xx[(min + 2) % DIM]) != 1) {
                    min = (min + 1) % DIM;
                    if (NOD(xx[(min + 1) % DIM], xx[(min + 2) % DIM]) != 1) {
                        error("error min element\n");
                    }
                }
            }

            gcd(xx[(min + 1) % DIM], xx[(min + 2) % DIM], u, v);
            if (abs(get_det_2_2(xx[(min + 1) % DIM], -v, xx[(min + 2) % DIM], u)) != 1) {
                error("error in Bezout values\n");
            }
            if (abs(0 * xx[min] + u * xx[(min + 1) % DIM] + v * xx[(min + 2) % DIM]) != 1) {
                error("inner product error\n");
            }

            two_col[min] = 0;
            two_col[(min + 1) % DIM] = u;
            two_col[(min + 2) % DIM] = v;
            if (two_col[0] == 0 && two_col[1] == 0) {
                for (int i = 0; i < DIM; i++) {
                    uni[i][0] = xx[i];
                }
                uni[0][1] = xx[2];
                uni[1][1] = 0;
                uni[2][1] = 0;
                uni[2][2] = 0;
                uni[0][2] = 0;
                uni[1][2] = 1;
                break;
            }
            else if (!((abs(two_col[0]) > 1 || abs(two_col[1]) > 1) && two_col[2] != 0)) {
                for (int i = 0; i < DIM; i++) {
                    uni[i][0] = xx[i];
                }
                uni[0][1] = -two_col[1];
                uni[1][1] = two_col[0];
                uni[2][1] = 0;
                uni[2][2] = 1;
                uni[0][2] = two_col[2] == 0 ? 0 : -u * v; 
                uni[1][2] = two_col[2] == 0 ? 0 : -u * v; 
                break;
            }
            else if (abs(two_col[0]) > 1 && two_col[2] != 0) {
                for (int i = 0; i < DIM; i++) {
                    uni[i][0] = xx[i];
                }
                uni[0][1] = -two_col[1];
                uni[1][1] = 1;
                uni[2][1] = 0;
                uni[2][2] = two_col[0];
                uni[0][2] = -two_col[2];
                uni[1][2] = 0;
                break;
            }
            else {
                for (int i = 0; i < DIM; i++) {
                    uni[i][0] = xx[i];
                }
                uni[0][1] = 1;
                uni[1][1] = two_col[0];
                uni[2][1] = 0;
                uni[2][2] = -two_col[1];
                uni[0][2] = 0;
                uni[1][2] = two_col[2];
                break;
            }

            if (abs(xx[min % DIM]) == abs(xx[(min + 1) % DIM])) {
                min = (min + 1) % DIM;
            }
            else if (abs(xx[min % DIM]) == abs(xx[(min + 2) % DIM])) {
                min = (min + 2) % DIM;
            }
            else
                min = abs(xx[(min + 1) % DIM]) < abs(xx[(min + 2) % DIM]) ? (min + 1) % DIM : (min + 2) % DIM;

        }

        // --_--
 
    }
    printf("multiplication by the first unimodular matrix\n");
    return uni;
}

std::vector<std::vector<int>> Algorithm::get_two_uni(const std::vector<std::vector<int>>& similar) {
    std::vector<std::vector<int>> submatrix_2_2(2, std::vector<int>(2));
    for (int i = 1; i < DIM; i++) {
        for (int j = 1; j < DIM; j++) {
            submatrix_2_2[i - 1][j - 1] = similar[i][j];
        }
    }

    int sub_eigenvalue_1 = 0, sub_eigenvalue_2 = 0, sub_eigenvalue = 0;
    int polynomial_A2 = -submatrix_2_2[0][0] - submatrix_2_2[1][1];
    int polynomial_A3 = submatrix_2_2[0][0] * submatrix_2_2[1][1] - submatrix_2_2[0][1] * submatrix_2_2[1][0];

    int determinant = sqrt((pow(polynomial_A2, 2) - 4 * polynomial_A3));
    if (pow(determinant, 2) == pow(polynomial_A2, 2) - 4 * polynomial_A3) {
        sub_eigenvalue_1 = (-polynomial_A2 - determinant) / 2;
        sub_eigenvalue_2 = (-polynomial_A2 + determinant) / 2;
    }
    else {
       error("not all integer eigenvalues\n");
    }

    sub_eigenvalue = std::min(sub_eigenvalue_1, sub_eigenvalue_2);
    std::vector<std::vector<int>> characteristic_matrix_A(2, std::vector<int>(2));
    characteristic_matrix_A = submatrix_2_2;
    for (int i = 0; i < 2; i++) {
        characteristic_matrix_A[i][i] -= sub_eigenvalue;
    }

    int nod;
    std::vector<int> eigenvector(2);
    if (characteristic_matrix_A[0][0] != 0 || characteristic_matrix_A[0][1] != 0) {
        nod = NOD(characteristic_matrix_A[0][0], characteristic_matrix_A[0][1]);
        eigenvector[0] = characteristic_matrix_A[0][0];
        eigenvector[1] = characteristic_matrix_A[0][1];
    }
    else if (characteristic_matrix_A[1][0] != 0 || characteristic_matrix_A[1][1] != 0) {
        nod = NOD(characteristic_matrix_A[1][0], characteristic_matrix_A[1][1]);
        eigenvector[0] = characteristic_matrix_A[1][0];
        eigenvector[1] = characteristic_matrix_A[1][1];
    }
    else {
        nod = 1;
        eigenvector[0] = 0;
        eigenvector[1] = 1;
    }
    while (nod != 1) {
        eigenvector[0] /= nod;
        eigenvector[1] /= nod;
        nod = NOD(eigenvector[0], eigenvector[1]);
    }
    std::vector<std::vector<int>> two_similar(DIM, std::vector<int>(DIM));
    two_similar[0][0] = 1;
    two_similar[1][1] = eigenvector[1];
    two_similar[2][1] = -eigenvector[0];
    int x, y;
    gcd(two_similar[1][1], two_similar[2][1], x, y);
    two_similar[1][2] = -y;
    two_similar[2][2] = x;
    printf("multiplication by the second unimodular matrix\n");
    return two_similar;
}

std::vector<int> Algorithm::gauss(const std::vector<std::vector<int>>& _matrix) {
    int rows = DIM;
    std::vector<std::vector<int>> gss;
    gss = _matrix;
    int k = 0;
    int d = 1;
    bool F = true;
    int iter = 0;
    for (int i = 0; i < DIM - iter; i++) {
        double k1 = gss[i][0] / (double)gss[(i + 1) % (DIM - iter)][0];
        double k2 = gss[i][1] / (double)gss[(i + 1) % (DIM - iter)][1];
        double k3 = gss[i][2] / (double)gss[(i + 1) % (DIM - iter)][2];
        if (i != 2) {
            if (gss[i][0] == 0 && gss[i][1] == 0 && gss[i][2] == 0) {
                for (int j = 0; j < DIM; j++) {
                    gss[i][j] = 0;
                    std::swap(gss[i][j], gss[DIM - 1 - iter][j]);
                }
                iter++;
                continue;
            }
            
        }

        if (k1 == k2 && k2 == k3) {
            for (int j = 0; j < DIM; j++) {
                gss[i][j] = 0;
                std::swap(gss[i][j], gss[2][j]);
            }
            iter++;
        }
    }
    for (int i = 0; i < DIM; i++) {
        int nd = (NOD(NOD(gss[i][0], gss[i][1]), gss[i][2]));
        for (int j = 0; j < DIM; j++) {
            if (nd != 0) gss[i][j] /= nd;
           // printf("%d  ", gss[i][j]);
        }
        //printf("\n");
    }

    while (k < DIM && F) {
        if (gss[k][k] == 0) {
            int i = k + 1;
            while (i < 3 && gss[i][k] == 0) {
                i++;
            }
            if (i < 3) {
                int j = k;
                while (j < 3) {
                    std::swap(gss[i][j], gss[k][j]);
                    j++;
                }
            }
            else {
                F = false;
                break;
            }
        }
        if (F) {
            int i = 0;
            while (i < DIM) {
                if (i != k) {
                    for (int j = k + 1; j < DIM; j++) {
                        int m = (gss[i][j] * gss[k][k] - gss[k][j] * gss[i][k]) % d;
                        if (m != 0) {
                            exit(-22);
                        }
                        gss[i][j] = (gss[i][j] * abs(gss[k][k]) - gss[k][j] * gss[i][k]) / d;
                    }
                    gss[i][k] = 0;
                    //gss[i][i] = gss[k][k];
                }
                i++;
            }
        }
        d = gss[k][k];
        if(gss[k + 1][k + 1] != 0)
        if (k + 1 < DIM - 1) gss[k][k] = abs(gss[k + 1][k + 1]);
        k++;

    }

    printf("gauss\n");
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            printf("%d  ", _matrix[i][j]);
        }
        printf("\t-->\t");
        int nd = (NOD(NOD(gss[i][0], gss[i][1]), gss[i][2]));
        for (int j = 0; j < DIM; j++) {
            if (nd != 0) gss[i][j] /= nd;
            printf("%d  ", gss[i][j]);
        }
        printf("\n");
    }
    bool flag = false;
    for (int i = 0; i < DIM; i++) {
        if (gss[i][i] == 0) {
            flag = true;
        }
    }
    if (!flag) {
        error("diagonal gauss\n");
    }
    std::vector<int> xx(DIM, -1);
    std::vector<int> x(DIM, -1);
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            if (gss[i][j % DIM] != 0 && gss[i][(j + 1) % DIM] == 0 && gss[i][(j + 2) % DIM] == 0) {
                x[j] = 0;
                break;
            }
        }
    }
    if (gss[1][0] == 0 && gss[1][1] == 0 && gss[1][2] == 0) {
        while (true) {
            for (int i = 0; i < DIM; i++) {
                if (/*x[i] == 0 && */ x[(i + 1) % DIM] != 0 && x[(i + 2) % DIM] != 0) {
                    xx[i] = 0;
                    xx[(i + 1) % DIM] = -gss[0][(i + 2) % DIM];
                    xx[(i + 2) % DIM] = gss[0][(i + 1) % DIM];
                    if (xx[(i + 1) % DIM] == 0 && xx[(i + 2) % DIM] == 0) {
                        xx[(i + 1) % DIM] = 1;
                        xx[(i + 2) % DIM] = 1;
                    }
                    goto down;
                }
            }
            error("can't find eigenvector\n");
        }
    }
    else if (gss[1][1] != 0 || gss[1][2] != 0) { // &&
        xx[2] = x[2] == 0 ? 0 : gss[1][1];
        /*if (-(gss[1][2] * xx[2]) % gss[1][1] != 0) {
            error("can't find eigenvector\n");
        }*/
        xx[1] = x[1] == 0 ? 0 : (xx[2] != 0 ? (-gss[1][2] * xx[2] / gss[1][1]) : gss[1][1] == 0 ? 1 : 0);
        if (gss[0][0] != 0) {
            if (-(gss[0][1] * xx[1] + gss[0][2] * xx[2]) % gss[0][0] != 0) {
                xx[0] = -(gss[0][1] * xx[1] + gss[0][2] * xx[2]);
                xx[1] *= abs(gss[0][0]);
                xx[2] *= abs(gss[0][0]);
            }
            else {
                xx[0] = -(gss[0][1] * xx[1] + gss[0][2] * xx[2]) / gss[0][0];
            }
        }
        else {
            xx[0] = 1;
        }
    }
    for (int i = 0; i < DIM; i++) {
        if (x[i] == -1 && x[(i + 1) % DIM] == 0 && x[(i + 2) % DIM] == 0) {
            xx[i] = 1;
            xx[(i + 1) % DIM] = 0;
            xx[(i + 2) % DIM] = 0;
        }
    }
down:
    int nod = 0;
    while (nod != 1) {
        nod = NOD(NOD(xx[0], xx[1]), xx[2]);
        xx[0] /= nod;
        xx[1] /= nod;
        xx[2] /= nod;
    }
    if (xx[0] * gss[0][0] + xx[1] * gss[0][1] + xx[2] * gss[0][2] != 0 || xx[0] * gss[1][0] + xx[1] * gss[1][1] + xx[2] * gss[1][2] != 0) {
        error("wrong eigenvector\n");
    }

    printf("eigenvector : (");
    printf("%d, ", xx[0]);
    printf("%d, ", xx[1]);
    printf("%d)^T\n\n", xx[2]);

    return xx;
}
