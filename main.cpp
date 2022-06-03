#include "Algorithm_matrix_dimension_3.h"

int main() {
	int input;
	printf("1 - find the canonical form\n2 - compare two matrices\n");
	scanf_s("%d", &input);
	system("cls");
	std::vector<std::vector<int>> A(DIM, std::vector<int>(DIM));
	std::vector<std::vector<int>> unimodular_A(DIM, std::vector<int>(DIM));
	for (int i = 0; i < DIM; i++) {
		unimodular_A[i][i] = 1;
	}
	if (input == 1) {
		/*A[0][0] = 2;
		A[0][1] = 2;
		A[0][2] = 3;
		A[1][0] = 2;
		A[1][1] = 5;
		A[1][2] = 6;
		A[2][0] = 1;
		A[2][1] = 1;
		A[2][2] = 6;*/

		/*A[0][0] = 2;
		A[0][1] = -21;
		A[0][2] = 10;
		A[1][0] = 0;
		A[1][1] = 8;
		A[1][2] = -17;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 8;*/

		/*A[0][0] = -1;
		A[0][1] = 4;
		A[0][2] = 232;
		A[1][0] = 0;
		A[1][1] = 1;
		A[1][2] = 0;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 1;*/

		/*A[0][0] = 5;
		A[0][1] = -1;
		A[0][2] = -1;
		A[1][0] = 0;
		A[1][1] = 4;
		A[1][2] = -1;
		A[2][0] = 0;
		A[2][1] = -1;
		A[2][2] = 4;*/

		/*A[0][0] = 5;
		A[0][1] = -1;
		A[0][2] = -1;
		A[1][0] = 0;
		A[1][1] = 5;
		A[1][2] = -1;
		A[2][0] = 0;
		A[2][1] = -1;
		A[2][2] = 5;*/

		/*A[0][0] = 5;
		A[0][1] = 6;
		A[0][2] = 3;
		A[1][0] = -1;
		A[1][1] = 0;
		A[1][2] = 1;
		A[2][0] = 1;
		A[2][1] = 2;
		A[2][2] = -1;*/

		/*A[0][0] = 1;
		A[0][1] = 5;
		A[0][2] = 6;
		A[1][0] = 0;
		A[1][1] = -2;
		A[1][2] = 1;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = -2;*/

		/*A[0][0] = 3;
		A[0][1] = 0;
		A[0][2] = -1;
		A[1][0] = 0;
		A[1][1] = 0;
		A[1][2] = 1;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 0;*/

		/*A[0][0] = 8;
		A[0][1] = 0;
		A[0][2] = 8;
		A[1][0] = 0;
		A[1][1] = 8;
		A[1][2] = 0;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 8;*/

		/*A[0][0] = 0;
		A[0][1] = 4;
		A[0][2] = 0;
		A[1][0] = 9;
		A[1][1] = 0;
		A[1][2] = 0;
		A[2][0] = 0;
		A[2][1] = 3;
		A[2][2] = 0;*/

		/*A[0][0] = 5;
		A[0][1] = 1;
		A[0][2] = 1;
		A[1][0] = 0;
		A[1][1] = 8;
		A[1][2] = 1;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 4;*/

		/*A[0][0] = 1;
		A[0][1] = 0;
		A[0][2] = 0;
		A[1][0] = 0;
		A[1][1] = 0;
		A[1][2] = 2;
		A[2][0] = 7;
		A[2][1] = 0;
		A[2][2] = 0;*/

		/*A[0][0] = 8;
		A[0][1] = 0;
		A[0][2] = 2;
		A[1][0] = 0;
		A[1][1] = -1;
		A[1][2] = 2;
		A[2][0] = 0;
		A[2][1] = 4;
		A[2][2] = 1;*/

		/*A[0][0] = 2;
		A[0][1] = 0;
		A[0][2] = 0;
		A[1][0] = 0;
		A[1][1] = 0;
		A[1][2] = 3;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 0;*/

		A[0][0] = 7;
		A[0][1] = 0;
		A[0][2] = 7;
		A[1][0] = 3;
		A[1][1] = 11;
		A[1][2] = 0;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 6;

		//A[0][0] = 1;
		//A[0][1] = 2;
		//A[0][2] = -7;
		//A[1][0] = 0;
		//A[1][1] = 1;
		//A[1][2] = 3;
		//A[2][0] = 0;
		//A[2][1] = 0;
		//A[2][2] = 1;

		/*A[0][0] = 5;
		A[0][1] = -3;
		A[0][2] = 0;
		A[1][0] = 0;
		A[1][1] = 3;
		A[1][2] = -2;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 0;*/

		/*A[0][0] = 6;
		A[0][1] = -14;
		A[0][2] = 0;
		A[1][0] = 0;
		A[1][1] = 21;
		A[1][2] = -6;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 0;*/

		//A[0][0] = 11;
		//A[0][1] = 12;
		//A[0][2] = 3;
		//A[1][0] = 0;
		//A[1][1] = 1;
		//A[1][2] = 7;
		//A[2][0] = 0;
		//A[2][1] = 0;
		//A[2][2] = 1;


		/*printf("input matrix :\n");
		for (int i = 0; i < DIM; i++) {
				scanf_s("%d %d %d", &A[i][0], &A[i][1], & A[i][2]);
			printf("\n");
		}
		system("cls");*/

		printf("input matrix :\n");
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++)
				A[i][j] < 0 ? printf("%d ", A[i][j]) : printf("%d  ", A[i][j]);
			printf("\n");
		}
		printf("\n");
		std::vector<std::vector<int>> canonical_matrix(DIM, std::vector<int>(DIM));
		Algorithm algo(A);
		canonical_matrix = algo.algorithm(A, unimodular_A);
		printf("matrix check\n");
		std::vector<std::vector<int>> check(DIM, std::vector<int>(DIM));
		std::vector<std::vector<int>> un({ { 1, 0, 0 }, {0, 1, 0}, {0, 0, 1} });

		check = multiplication_and_inverse_B(A, unimodular_A, un);
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++) {
				if (canonical_matrix[i][j] != check[i][j]) {
					printf("check failed\n");
					return -404;
				}
			}
		}
		printf("check done\n");
		system("pause");
	}
	else {
		std::vector<std::vector<int>> B(DIM, std::vector<int>(DIM));
		std::vector<std::vector<int>> unimodular_B(DIM, std::vector<int>(DIM));
		for (int i = 0; i < DIM; i++) {
			unimodular_B[i][i] = 1;
		}

		/*printf("input matrix A:\n");
		for (int i = 0; i < DIM; i++) {
				scanf_s("%d %d %d", &A[i][0], &A[i][1], &A[i][2]);
			printf("\n");
		}
		system("cls");

		printf("input matrix B:\n");
		for (int i = 0; i < DIM; i++) {
				scanf_s("%d %d %d", &B[i][0], &B[i][1], &B[i][2]);
			printf("\n");
		}
		system("cls");*/

		A[0][0] = 6;
		A[0][1] = -14;
		A[0][2] = 0;
		A[1][0] = 0;
		A[1][1] = 21;
		A[1][2] = -6;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 0;

		B[0][0] = 0;
		B[0][1] = 2;
		B[0][2] = 20;
		B[1][0] = 0;
		B[1][1] = 6;
		B[1][2] = 6;
		B[2][0] = 0;
		B[2][1] = 0;
		B[2][2] = 21;


		std::vector<std::vector<int>> canonical_matrix_A(DIM, std::vector<int>(DIM));
		Algorithm algo_a(A);
		canonical_matrix_A = algo_a.algorithm(A, unimodular_A);
		printf("matrix check\n");
		std::vector<std::vector<int>> check_A(DIM, std::vector<int>(DIM));
		std::vector<std::vector<int>> un_A({ { 1, 0, 0 }, {0, 1, 0}, {0, 0, 1} });

		check_A = multiplication_and_inverse_B(A, unimodular_A, un_A);
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++) {
				if (canonical_matrix_A[i][j] != check_A[i][j]) {
					printf("check failed\n");
					return -404;
				}
			}
		}

		std::vector<std::vector<int>> canonical_matrix_B(DIM, std::vector<int>(DIM));
		Algorithm algo_b(B);
		canonical_matrix_B = algo_b.algorithm(B, unimodular_B);
		printf("matrix check\n");
		std::vector<std::vector<int>> check_B(DIM, std::vector<int>(DIM));
		std::vector<std::vector<int>> un_B({ { 1, 0, 0 }, {0, 1, 0}, {0, 0, 1} });

		check_B = multiplication_and_inverse_B(B, unimodular_B, un_B);
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++) {
				if (canonical_matrix_B[i][j] != check_B[i][j]) {
					printf("check failed\n");
					return -404;
				}
			}
		}

		system("cls");

		printf("input matrix A:\n");
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++)
				A[i][j] < 0 ? printf("%d ", A[i][j]) : printf("%d  ", A[i][j]);
			printf("\n");
		}
		printf("\n");

		printf("input matrix B:\n");
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++)
				B[i][j] < 0 ? printf("%d ", B[i][j]) : printf("%d  ", B[i][j]);
			printf("\n");
		}
		printf("\n");

		printf("matrix check\n");
		std::vector<std::vector<int>> check_AB(DIM, std::vector<int>(DIM));
		std::vector<std::vector<int>> un_AB({ { 1, 0, 0 }, {0, 1, 0}, {0, 0, 1} });

		check_AB = multiplication_and_inverse_B(B, multiplication(unimodular_B, inverse_matrix(unimodular_A)), un_AB);
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++) {
				if (A[i][j] != check_AB[i][j]) {
					printf("false\n");
					return -404;
				}
			}
		}

		printf("true\n\n");
		system("pause");
	}
	return 0;
}