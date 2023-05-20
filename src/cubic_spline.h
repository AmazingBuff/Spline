//
// Created by Amazing on 2023/5/19.
//

#pragma once

#include"matrix_decompose.h"
#include<cmath>
#include<cassert>


namespace Amazing
{
    //basic function
    //f(x) = a + b(x - x0) + c(x - x0)^2 + d(x - x0)^3
    //f'(x) = b + 2c(x - x0) + 3d(x - x0)^2
    //f''(x) = 2c + 6(x - x0)

	struct Coefficient
	{
		std::vector<double> a;
		std::vector<double> b;
		std::vector<double> c;
		std::vector<double> d;

		explicit Coefficient(const size_t size)
		{
			a.resize(size);
			b.resize(size);
			c.resize(size);
			d.resize(size);
		}
	};


    class CubicSpline
    {
    public:
		enum class BoundaryConditionType
		{
			Natural,
			Clamped,
			Not_A_Knot
		};

		[[nodiscard]] static Coefficient interpolation(const std::vector<double>& x, const std::vector<double>& y,
		BoundaryConditionType type = BoundaryConditionType::Natural, double a = 0.0f, double b = 0.0f)
		{
			assert(x.size() == y.size());

			static size_t dist = x.size();

			Matrix A(dist);

			std::vector<double> B(dist);
			memset(B.data(), 0, sizeof(double) * dist);

			for(size_t i = 1; i < dist - 1; i++)
			{
				double h_i_1 = x[i] - x[i - 1];
				double h_i = x[i + 1] - x[i];
				A[i][i - 1] = h_i_1;
				A[i][i] = 2 * (h_i_1 + h_i);
				A[i][i + 1] = h_i;

				B[i] = (y[i + 1] - y[i]) / h_i - (y[i] - y[i - 1]) / h_i_1;
			}

			switch (type)
			{
				case BoundaryConditionType::Natural:
				{
					//natural boundary condition
					//S_0''(x_0) = S_n_1''(x_n) = 0 -->m_0 = 0, m_n = 0
					A[0][0] = 1.0f;
					A[dist - 1][dist - 1] = 1.0f;

					break;
				}
				case BoundaryConditionType::Clamped:
				{
					//clamped boundary condition
					//S_0'(x_0) = a, S_n_1'(x_n) = b
					double h_0 = x[1] - x[0];
					double h_n_1 = x[dist - 1] - x[dist - 2];
					A[0][0] = 2 * h_0;
					A[0][1] = h_0;
					B[0] = 6 * ((y[1] - y[0]) / h_0 - a);
					A[dist - 1][dist - 2] = h_n_1;
					A[dist - 1][dist - 1] = 2 * h_n_1;
					B[dist - 1] = 6 * (b - (y[dist - 1] - y[dist - 2]) / h_n_1);

					break;
				}
				case BoundaryConditionType::Not_A_Knot:
				{
					//not a knot
					//S_0'''(x_0) = S_1'''(x_1), S_n_2'''(x_n_1) = S_n_1'''(x_n)
					double h_0 = x[1] - x[0];
					double h_1 = x[2] - x[1];
					double h_n_1 = x[dist - 1] - x[dist - 2];
					double h_n_2 = x[dist - 2] - x[dist - 3];
					A[0][0] = -h_1;
					A[0][1] = h_0 + h_1;
					A[0][2] = -h_0;
					A[dist - 1][dist - 3] = -h_n_1;
					A[dist - 1][dist - 2] = -h_n_1 - h_n_2;
					A[dist - 1][dist - 3] = -h_n_2;

					break;
				}
			}

			TriangularMatrix l_matrix(dist);
			TriangularMatrix u_matrix(dist);
			SquareMatrixDecompose::LUDecompose(A, l_matrix, u_matrix);

			std::vector<double> m = u_matrix.inverse_upper() * (l_matrix.inverse_lower() * B);

			Coefficient ret(dist - 1);
			for(size_t i = 0; i < dist - 1; i++)
			{
				double h_i = x[i + 1] - x[i];
				double dm_i = m[i + 1] - m[i];

				ret.a[i] = y[i];
				ret.b[i] = (y[i + 1] - y[i]) / h_i - h_i * m[i] / 2.0f - h_i * dm_i / 6.0f;
				ret.c[i] = m[i] / 2.0f;
				ret.d[i] = dm_i / (6.0f * h_i);
			}

			return ret;
		}

//        void interpolation()
//        {
//            size_t dist = x.size() - 1;
//
//            //[0, dist) ai
//            //[dist, 2 * dist) bi
//            //[2 * dist, 3 * dist) ci
//            //[3 * dist, 4 * dist) di
//            std::vector<std::vector<double>> A;
//            A.resize(dist * 4);
//            for(size_t i = 0; i < A.size(); i++)
//            {
//                A[i].resize(dist * 4);
//                memset(A[i].data(), 0, sizeof(double) * dist * 4);
//            }
//
//            std::vector<double> Y(dist * 4);
//            memset(Y.data(), 0, sizeof(double) * dist * 4);
//            //Si(xi) = yi
//            for(size_t i = 0; i < dist; i++)
//            {
//                A[i][i] = 1.0f;
//                Y[i] = y[i];
//            }
//
//            //Si(xi+1) = yi+1
//            for(size_t i = 0; i < dist; i++)
//            {
//                double hi = x[i + 1] - x[i];
//
//                A[i + dist][4 * i] = 1.0f;
//                A[i + dist][i + dist] = hi;
//                A[i + dist][i + 2 * dist] = std::pow(hi, 2);
//                A[i + dist][i + 3 * dist] = std::pow(hi, 3);
//                Y[i + dist] = y[i + 1];
//            }
//
//            //Si'(xi+1) = Si+1'(xi+1)
//            for(size_t i = 0; i < dist - 1; i++)
//            {
//                double hi = x[i + 1] - x[i];
//
//                A[i + 2 * dist][i + dist] = 1.0f;
//                A[i + 2 * dist][i + 1 + dist] = -1.0f;
//                A[i + 2 * dist][i + 2 * dist] = 2.0f * hi;
//                A[i + 2 * dist][i + 3 * dist] = 3.0f * std::pow(hi, 2);
//            }
//
//            //Si''(xi+1) = Si+1''(xi+1)
//            for(size_t i = 0; i < dist - 1; i++)
//            {
//                double hi = x[i + 1] - x[i];
//
//                A[i + 3 * dist][i + 2 * dist] = 2.0f;
//                A[i + 3 * dist][i + 1 + 2 * dist] = -2.0f;
//                A[i + 3 * dist][i + 3 * dist] = 6.0f * hi;
//            }
//
//            //natural boundary condition
//            //S0''(x0) = 0
//            A[4 * dist - 2][2 * dist] = 2.0f;
//            //Sn''(xn) = 0
//            A[4 * dist - 1][3 * dist - 1] = 2.0f;
//        }

    private:
        CubicSpline() = default;
    };


    class EquidistanceCubicSpline
    {
        EquidistanceCubicSpline();
    };
}