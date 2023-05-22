//
// Created by Amazing on 2023/5/22.
//

#pragma once

#include<cstddef>
#include<initializer_list>
#include<algorithm>

namespace Amazing
{
    class Vector
    {
    public:
        Vector() : m_size(0), m_vector(nullptr) {}

        ~Vector()
        {
            if(m_vector)
                clear();
        }

        explicit Vector(const size_t m) : m_size(m)
        {
            m_vector = new double[m];
            memset(m_vector, 0, sizeof(double) * m);
        }

        Vector(const double* begin, const double* end)
        {
            m_size = end - begin;
            m_vector = new double[m_size];
            memcpy(m_vector, begin, sizeof(double) * m_size);
        }

        Vector(const std::initializer_list<double>& list)
        {
            m_size = list.size();
            m_vector = new double[m_size];
            for(auto i : list)
                *m_vector++ = i;
        }

        Vector(const Vector& vector)
        {
            m_size = vector.size();
            m_vector = new double[m_size];
            memcpy(m_vector, &vector[0], sizeof(double) * m_size);
        }

        Vector& operator=(const Vector& other)
        {
            if(this != &other)
            {
                delete[] m_vector;
                m_size = other.size();
                m_vector = new double[m_size];
                memcpy(m_vector, other.data(), sizeof(double) * m_size);
            }
            return *this;
        }

        double& operator[](size_t index)
        {
            return m_vector[index];
        }

        const double& operator[](size_t index) const
        {
            return m_vector[index];
        }

        double operator*(const Vector& other) const
        {
            double ret = 0.0f;
            for(size_t i = 0; i < std::min(m_size, other.size()); i++)
                ret += m_vector[i] * other[i];
            return ret;
        }

        [[nodiscard]] size_t size() const
        {
            return m_size;
        }

        [[nodiscard]] double* data() const
        {
            return m_vector;
        }

        void clear()
        {
            delete[] m_vector;
            m_vector = nullptr;
            m_size = 0;
        }

        void reset(const size_t m)
        {
            if(m != m_size)
            {
                clear();
                m_vector = new double[m];
                memset(m_vector, 0, sizeof(double) * m);
                m_size = m;
            }
        }

    private:
        size_t m_size;
        double* m_vector;
    };
}
