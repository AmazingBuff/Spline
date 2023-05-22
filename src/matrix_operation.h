//
// Created by AmazingBuff on 2023/5/20.
//

#pragma once

#include"vector_operation.h"
#include<array>
#include<cstddef>
#include<cstring>

namespace Amazing
{
	class Matrix
	{
	public:
        Matrix() : m_row(0), m_col(0), m_matrix(nullptr) {}

        ~Matrix()
        {
            if(m_matrix)
                clear();
        }

		Matrix(const size_t m, const size_t n) : m_row(m), m_col(n)
		{
			m_matrix = new Vector[m];
			for(size_t i = 0; i < m; i++)
				m_matrix[i].reset(n);
		}

		explicit Matrix(const size_t m) : m_row(m), m_col(m)
		{
            m_matrix = new Vector[m];
            for(size_t i = 0; i < m; i++)
                m_matrix[i].reset(m);
		}

        Matrix(const Matrix& matrix)
        {
            m_row = matrix.row();
            m_col = matrix.col();
            m_matrix = new Vector[m_row];
            for(size_t i = 0; i < m_row; i++)
            {
                m_matrix[i].reset(m_col);
                memcpy(m_matrix[i].data(), matrix[i].data(), sizeof(double) * m_col);
            }
        }

        Matrix(const std::initializer_list<std::initializer_list<double>>& lists)
        {
            m_row = lists.size();
            m_matrix = new Vector[m_row];

            size_t i = 0;
            for(auto& list : lists)
            {
                m_matrix[i] = Vector(list.begin(), list.end());
                i++;
            }


            m_col = m_matrix[0].size();
        }

		const Vector& operator[](size_t index) const
		{
			return m_matrix[index];
		}

        Vector& operator[](size_t index)
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

		[[nodiscard]] Vector operator*(const Vector& vector) const
		{
			Vector ret(m_row);

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

		[[nodiscard]] size_t row() const
		{
			return m_row;
		}

        [[nodiscard]] size_t col() const
		{
			return m_col;
		}

        void clear()
        {
            for(size_t i = 0; i < m_row; i++)
                m_matrix[i].clear();
            delete[] m_matrix;
            m_matrix = nullptr;

            m_row = 0;
            m_col = 0;
        }

        void reset(const size_t m, const size_t n)
        {
            if(m != m_row || n != m_col)
            {
                clear();
                m_matrix = new Vector[m];
                for(size_t i = 0; i < m; i++)
                    m_matrix[i].reset(n);

                m_row = m;
                m_col = n;
            }
        }

        void reset(const size_t m)
        {
            if(m != m_row || m != m_col)
            {
                clear();
                m_matrix = new Vector[m];
                for(size_t i = 0; i < m; i++)
                    m_matrix[i].reset(m);
                m_row = m;
                m_col = m;
            }
        }

	protected:
		size_t m_row;
		size_t m_col;
		Vector* m_matrix;
	};


	class TriangularMatrix : public Matrix
	{
	public:
		explicit TriangularMatrix(const size_t m) : Matrix(m) {}

        TriangularMatrix(const std::initializer_list<std::initializer_list<double>>& list) : Matrix(list) {}

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
