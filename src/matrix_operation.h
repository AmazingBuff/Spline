//
// Created by AmazingBuff on 2023/5/20.
//

#pragma once

#include<vector>
#include<array>
#include<cstddef>
#include<cstring>

namespace Amazing
{
	class Matrix
	{
	public:
		Matrix(const size_t m, const size_t n) : m_row(m), m_col(n)
		{
			m_matrix = new double*[m];
			for(size_t i = 0; i < m; i++)
			{
				m_matrix[i] = new double[n];
				memset(&m_matrix[i][0], 0, sizeof(double) * n);
			}
		}

		explicit Matrix(const size_t m) : m_row(m), m_col(m)
		{
			m_matrix = new double*[m];
			for(size_t i = 0; i < m; i++)
			{
				m_matrix[i] = new double[m];
				memset(&m_matrix[i][0], 0, sizeof(double) * m);
			}
		}

		explicit Matrix(const std::vector<std::vector<double>>& matrix)
		{
			m_row = matrix.size();
			m_col = matrix[0].size();
			m_matrix = new double*[m_row];
			for(size_t i = 0; i < m_row; i++)
			{
				m_matrix[i] = new double[m_col];
				memcpy(&m_matrix[i][0], matrix[i].data(), sizeof(double) * m_col);
			}
		}

		const double* operator[](size_t index) const
		{
			return m_matrix[index];
		}

		double* operator[](size_t index)
		{
			return m_matrix[index];
		}

		void transpose()
		{
			for(size_t i = 0; i < m_row; i++)
			{
				for(size_t j = i + 1; j < m_row; j++)
				{
					double tmp = m_matrix[i][j];
					m_matrix[i][j] = m_matrix[j][i];
					m_matrix[j][i] = tmp;
				}
			}
		}

		[[nodiscard]] std::vector<double> operator*(const std::vector<double>& vector) const
		{
			std::vector<double> ret(m_row);
			memset(ret.data(), 0, sizeof(double) * m_row);

			size_t intern = std::min(m_row, vector.size());
			for(size_t i = 0; i < m_row; i++)
			{
				for(size_t j = 0; j < intern; j++)
					ret[i] += m_matrix[i][j] * vector[j];
			}

			return ret;
		}

		[[nodiscard]] Matrix operator*(const Matrix& matrix) const
		{
			Matrix ret(m_row);
			for(size_t i = 0; i < m_row; i++)
			{
				for(size_t j = 0; j < m_row; j++)
				{
					for(size_t k = 0; k < m_row; k++)
						ret[i][j] += m_matrix[i][k] * matrix[k][j];
				}
			}
			return ret;
		}

		size_t row() const
		{
			return m_row;
		}

		size_t col() const
		{
			return m_col;
		}

	protected:
		size_t m_row;
		size_t m_col;
		double** m_matrix;
	};


	class TriangularMatrix : public Matrix
	{
	public:
		explicit TriangularMatrix(const size_t m) : Matrix(m) {}

		explicit TriangularMatrix(const std::vector<std::vector<double>>& matrix) : Matrix(matrix) {}

		[[nodiscard]] TriangularMatrix inverse_upper() const
		{
			TriangularMatrix ret(m_row);
			for(size_t i = 0; i < m_row; i++)
			{
				ret[i][i] = 1.0f / this->m_matrix[i][i];
				for(size_t j = i - 1; j < i; j--)
				{
					double sum = 0.0f;
					for(size_t k = j + 1; k <= i; k++)
						sum += this->m_matrix[j][k] * ret[k][i];
					ret[j][i] = -sum * ret[j][j];
				}
			}

			return ret;
		}

		[[nodiscard]] TriangularMatrix inverse_lower() const
		{
			TriangularMatrix ret(m_row);
			for(size_t i = m_row - 1; i < m_row; i--)
			{
				ret[i][i] = 1.0f / this->m_matrix[i][i];
				for(size_t j = i + 1; j < m_row; j++)
				{
					double sum = 0.0f;
					for(size_t k = i; k < j; k++)
						sum += this->m_matrix[j][k] * ret[k][i];
					ret[j][i] = -sum * ret[j][j];
				}
			}

			return ret;
		}
	};
}
