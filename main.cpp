//
// Created by AmazingBuff on 2023/5/20.
//
#include<cubic_spline.h>

int main()
{
	std::vector<std::vector<double>> array =
		{
			{4,2,1,5},
			{8,7,2,10},
			{4,8,3,6},

			{6,8,4,9}
		};

	Amazing::Matrix matrix(array);
	Amazing::TriangularMatrix l_matrix(4);
	Amazing::TriangularMatrix u_matrix(4);
	Amazing::SquareMatrixDecompose::LUDecompose(matrix, l_matrix, u_matrix);


	std::vector<std::vector<double>> maa =
		{
			{2,2,-2},
			{0,3,-6},
			{0,0,-5},
		};
	Amazing::TriangularMatrix ma(maa);
	Amazing::TriangularMatrix upper = ma.inverse_upper();
	ma.transpose();
	Amazing::TriangularMatrix lower = ma.inverse_lower();

	std::vector<double> x = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
	std::vector<double> y = {1.1f, 2.6f, 3.9f, 7.4f, 9.3f, 13.6f};

	Amazing::Coefficient coef = Amazing::CubicSpline::interpolation(x, y);

	double f_x = coef.a[1] + coef.b[1] * (2.5f - 2.0f) + coef.c[1] * std::pow((2.5f - 2.0f), 2.0f) + coef.d[1] * std::pow((2.5f - 2.0f), 3.0f);

	system("pause");
}