//
// Created by Amazing on 2023/5/19.
//

#pragma once

#include<vector>
#include<cmath>


namespace Amazing
{
    //basic function
    //f(x) = a + b(x - x0) + c(x - x0)^2 + d(x - x0)^3
    //f'(x) = b + 2c(x - x0) + 3d(x - x0)^2
    //f''(x) = 2c + 6(x - x0)

    class CubicSpline
    {
    public:
        CubicSpline(const std::vector<double>& x, const std::vector<double>& y) : x(x), y(y)
        {
            assert(x.size() == y.size());
        }

        void interpolation()
        {
            size_t dist = x.size() - 1;

            //[0, dist) ai
            //[dist, 2 * dist) bi
            //[2 * dist, 3 * dist) ci
            //[3 * dist, 4 * dist) di
            std::vector<std::vector<double>> A;
            A.resize(dist * 4);
            for(size_t i = 0; i < A.size(); i++)
            {
                A[i].resize(dist * 4);
                memset(A[i].data(), 0, sizeof(double) * dist * 4);
            }

            std::vector<double> Y(dist * 4);
            memset(Y.data(), 0, sizeof(double) * dist * 4);
            //Si(xi) = yi
            for(size_t i = 0; i < dist; i++)
            {
                A[i][i] = 1.0f;
                Y[i] = y[i];
            }

            //Si(xi+1) = yi+1
            for(size_t i = 0; i < dist; i++)
            {
                double hi = x[i + 1] - x[i];

                A[i + dist][4 * i] = 1.0f;
                A[i + dist][i + dist] = hi;
                A[i + dist][i + 2 * dist] = std::pow(hi, 2);
                A[i + dist][i + 3 * dist] = std::pow(hi, 3);
                Y[i + dist] = y[i + 1];
            }

            //Si'(xi+1) = Si+1'(xi+1)
            for(size_t i = 0; i < dist - 1; i++)
            {
                double hi = x[i + 1] - x[i];

                A[i + 2 * dist][i + dist] = 1.0f;
                A[i + 2 * dist][i + 1 + dist] = -1.0f;
                A[i + 2 * dist][i + 2 * dist] = 2.0f * hi;
                A[i + 2 * dist][i + 3 * dist] = 3.0f * std::pow(hi, 2);
            }

            //Si''(xi+1) = Si+1''(xi+1)
            for(size_t i = 0; i < dist - 1; i++)
            {
                double hi = x[i + 1] - x[i];

                A[i + 3 * dist][i + 2 * dist] = 2.0f;
                A[i + 3 * dist][i + 1 + 2 * dist] = -2.0f;
                A[i + 3 * dist][i + 3 * dist] = 6.0f * hi;
            }

            //natural boundary condition
            //S0''(x0) = 0
            A[4 * dist - 2][2 * dist] = 2.0f;
            //Sn''(xn) = 0
            A[4 * dist - 1][3 * dist - 1] = 2.0f;
        }

    private:
        const std::vector<double>& x;
        const std::vector<double>& y;
    };


    class EquidistanceCubicSpline
    {
        EquidistanceCubicSpline();
    };
}