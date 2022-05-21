#include "Algorithm_matrix_dimension_3.h"

int main() {
	std::vector<std::vector<int>> vector(DIM, std::vector<int>(DIM));
	std::vector<std::vector<int>> analog(DIM, std::vector<int>(DIM));
	for (int i = 0; i < DIM; i++) {
		analog[i][i] = 1;
	}
	/*vector[0][0] = 2;
	vector[0][1] = 2;
	vector[0][2] = 3;
	vector[1][0] = 2;
	vector[1][1] = 5;
	vector[1][2] = 6;
	vector[2][0] = 1;
	vector[2][1] = 1;
	vector[2][2] = 6;*/

	/*vector[0][0] = 2;
	vector[0][1] = -21;
	vector[0][2] = 10;
	vector[1][0] = 0;
	vector[1][1] = 8;
	vector[1][2] = -17;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 8;*/

	/*vector[0][0] = -1;
	vector[0][1] = 4;
	vector[0][2] = 232;
	vector[1][0] = 0;
	vector[1][1] = 1;
	vector[1][2] = 0;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 1;*/

	/*vector[0][0] = 5;
	vector[0][1] = -1;
	vector[0][2] = -1;
	vector[1][0] = 0;
	vector[1][1] = 4;
	vector[1][2] = -1;
	vector[2][0] = 0;
	vector[2][1] = -1;
	vector[2][2] = 4;*/

	/*vector[0][0] = 5;
	vector[0][1] = -1;
	vector[0][2] = -1;
	vector[1][0] = 0;
	vector[1][1] = 5;
	vector[1][2] = -1;
	vector[2][0] = 0;
	vector[2][1] = -1;
	vector[2][2] = 5;*/

	/*vector[0][0] = 5;
	vector[0][1] = 6;
	vector[0][2] = 3;
	vector[1][0] = -1;
	vector[1][1] = 0;
	vector[1][2] = 1;
	vector[2][0] = 1;
	vector[2][1] = 2;
	vector[2][2] = -1;*/

	/*vector[0][0] = 1;
	vector[0][1] = 5;
	vector[0][2] = 6;
	vector[1][0] = 0;
	vector[1][1] = -2;
	vector[1][2] = 1;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = -2;*/

	/*vector[0][0] = 3;
	vector[0][1] = 0;
	vector[0][2] = -1;
	vector[1][0] = 0;
	vector[1][1] = 0;
	vector[1][2] = 1;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 0;*/

	/*vector[0][0] = 8;
	vector[0][1] = 0;
	vector[0][2] = 8;
	vector[1][0] = 0;
	vector[1][1] = 8;
	vector[1][2] = 0;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 8;*/

	/*vector[0][0] = 0;
	vector[0][1] = 4;
	vector[0][2] = 0;
	vector[1][0] = 9;
	vector[1][1] = 0;
	vector[1][2] = 0;
	vector[2][0] = 0;
	vector[2][1] = 3;
	vector[2][2] = 0;*/

	/*vector[0][0] = 5;
	vector[0][1] = 1;
	vector[0][2] = 1;
	vector[1][0] = 0;
	vector[1][1] = 8;
	vector[1][2] = 1;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 4;*/

	/*vector[0][0] = 1;
	vector[0][1] = 0;
	vector[0][2] = 0;
	vector[1][0] = 0;
	vector[1][1] = 0;
	vector[1][2] = 2;
	vector[2][0] = 7;
	vector[2][1] = 0;
	vector[2][2] = 0;*/

	/*vector[0][0] = 8;
	vector[0][1] = 0;
	vector[0][2] = 2;
	vector[1][0] = 0;
	vector[1][1] = -1;
	vector[1][2] = 2;
	vector[2][0] = 0;
	vector[2][1] = 4;
	vector[2][2] = 1;*/

	/*vector[0][0] = 2;
	vector[0][1] = 0;
	vector[0][2] = 0;
	vector[1][0] = 0;
	vector[1][1] = 0;
	vector[1][2] = 3;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 0;*/

	vector[0][0] = 7;
	vector[0][1] = 0;
	vector[0][2] = 7;
	vector[1][0] = 3;
	vector[1][1] = 11;
	vector[1][2] = 0;
	vector[2][0] = 0;
	vector[2][1] = 0;
	vector[2][2] = 6;

	//vector[0][0] = 1;
	//vector[0][1] = 2;
	//vector[0][2] = -7;
	//vector[1][0] = 0;
	//vector[1][1] = 1;
	//vector[1][2] = 3;
	//vector[2][0] = 0;
	//vector[2][1] = 0;
	//vector[2][2] = 1;


	/*printf("input matrix :\n");
	for (int i = 0; i < DIM; i++) {
			scanf_s("%d %d %d", &vector[i][0], &vector[i][1], & vector[i][2]);
		printf("\n");
	}
	system("cls");*/

	printf("input matrix :\n");
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++)
			vector[i][j] < 0 ? printf("%d ", vector[i][j]) : printf("%d  ", vector[i][j]);
		printf("\n");
	}
	printf("\n");
	std::vector<std::vector<int>> canonical_matrix(DIM, std::vector<int>(DIM));
	Algorithm algo(vector);
	canonical_matrix = algo.algorithm(vector, analog);
	printf("matrix check\n");
	std::vector<std::vector<int>> check(DIM, std::vector<int>(DIM));
	std::vector<std::vector<int>> un(DIM, std::vector<int>(DIM));
	for (int i = 0; i < DIM; i++) {
		un[i][i] = 1;
	}
	check = algo.multiplication_and_inverse_B(vector, analog, un);
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
	return 0;
}
