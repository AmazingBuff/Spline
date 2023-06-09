//
// Created by AmazingBuff on 2023/5/19.
//

#pragma once

#include"matrix_operation.h"
#include<chrono>
#include<iostream>

namespace Amazing
{
	class SquareMatrixDecompose
	{
	public:
        //diagonal element can't be zero
        //matrix, l_matrix and u_matrix: m x m
		static void LUDecompose(const Matrix& matrix, TriangularMatrix& l_matrix, TriangularMatrix& u_matrix)
		{
			size_t dist = matrix.row();
            if(l_matrix.col() != dist || l_matrix.row() != dist)
                l_matrix.reset(dist);
            if(u_matrix.col() != dist || u_matrix.row() != dist)
                u_matrix.reset(dist);

			//LU decompose
			for(size_t i = 0; i < dist; i++)
			{
				for(size_t j = i; j < dist; j++)
				{
					double u_post = 0.0f;
					for(size_t k = 0; k < i; k++)
						u_post += l_matrix[i][k] * u_matrix[k][j];
					u_matrix[i][j] = matrix[i][j] - u_post;
				}

				l_matrix[i][i] = 1.0f;
				for(size_t j = i + 1; j < dist; j++)
				{
					double l_post = 0.0f;
					for(size_t k = 0; k < i; k++)
						l_post += l_matrix[j][k] * u_matrix[k][i];
					l_matrix[j][i] = (matrix[j][i] - l_post) / u_matrix[i][i];
				}
			}


			for(size_t i = 0; i < dist; i++)
			{
				for(size_t j = 0; j < dist; j++)
					std::cout << l_matrix[i][j] << ' ';
				std::cout << std::endl;
			}

			for(size_t i = 0; i < dist; i++)
			{
				for(size_t j = 0; j < dist; j++)
					std::cout << u_matrix[i][j] << ' ';
				std::cout << std::endl;
			}
		}
	private:
		SquareMatrixDecompose() = default;
	};


    class MatrixDecompose
    {
    public:
        //matrix: m x n
        static void QRDecompose(const Matrix& matrix, Matrix& q, Matrix& r)
        {

        }

    private:
        MatrixDecompose() = default;
    };
}

