#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#define DIM 3
class Algorithm{
	std::vector<std::vector<int>>matrix;
	int first_coefficient,
		second_coefficient,
		third_coefficient,
		fourth_coefficient;
	int first_eigenvalues,
		second_eigenvalues,
		third_eigenvalues;
public:
	Algorithm(const std::vector<std::vector<int>>& matr) { matrix = matr; }
	std::vector<std::vector<int>> algorithm(const std::vector<std::vector<int>>& matr, std::vector<std::vector<int>>& analog);
	void get_eigenvalues(const std::vector<std::vector<int>>& matr);
	void print_coeffs();
	std::vector<std::vector<int>> get_unimod_matrix(const std::vector<std::vector<int>>& uni);
	std::vector<std::vector<int>> get_two_uni(const std::vector < std::vector<int>>& similar);
	std::vector<int> gauss(const std::vector<std::vector<int>>& matrix);
	int get_multiplicity(const int elem);
	void error(std::string str) {
		//system("cls");
		printf("%s", str.c_str());
		system("pause");
		exit(404);
	}
};

int NOD(int a, int b);
int gcd(int a, int b, int& x, int& y);
int get_det_2_2(const int& a, const int& b, const int& c, const int& d);
int get_det_3_3(const std::vector<std::vector<int>>& matr);
std::vector<std::vector<int>> multiplication(const std::vector<std::vector<int>>& a, const std::vector<std::vector<int>>& b);
std::vector<std::vector<int>> multiplication_and_inverse_B(const std::vector<std::vector<int>>& a, const std::vector<std::vector<int>>& b, std::vector<std::vector<int>>& analog);
std::vector<std::vector<int>> inverse_matrix(const std::vector<std::vector<int>>& matr);
bool pair(int a, int b, int c);
bool no_pair(int a, int b, int c);
int get_integer(int a, int b);