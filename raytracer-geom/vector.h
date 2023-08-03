#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <algorithm>
#include <cassert>

const double kEps = 0.00000001;

class Vector {
public:
    Vector() {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] = 0;
        }
    }
    Vector(const std::initializer_list<double>& list) {
        size_t ind = 0;
        for (const auto& elem : list) {
            data_[ind++] = elem;
        }
    }
    Vector(const std::array<double, 3>& data) : data_(data) {
    }

    double& operator[](size_t ind) {
        assert(ind < 3);
        return data_[ind];
    }
    double operator[](size_t ind) const {
        assert(ind < 3);
        return data_[ind];
    }

    void Normalize() {
        double tmp = 0;
        for (size_t i = 0; i < 3; ++i) {
            tmp += (data_[i] * data_[i]);
        }
        tmp = sqrt(tmp);
        for (size_t i = 0; i < 3; ++i) {
            data_[i] /= tmp;
        }
    }

    Vector& operator+=(const Vector& other) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    Vector operator+(const Vector& other) const {
        Vector ans = *this;
        ans += other;
        return ans;
    }

    Vector& operator-=(const Vector& other) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    Vector operator-(const Vector& other) const {
        Vector ans = *this;
        ans -= other;
        return ans;
    }

    Vector& operator*=(double scalar) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    Vector operator*(double scalar) const {
        Vector ans = *this;
        ans *= scalar;
        return ans;
    }

    //  Multiplication by elements
    Vector& operator*=(const Vector& other) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] *= other.data_[i];
        }
        return *this;
    }

    Vector operator*(const Vector& other) const {
        Vector ans = *this;
        ans *= other;
        return ans;
    }

    Vector& operator/=(double scalar) {
        assert(std::abs(scalar) > kEps);
        for (size_t i = 0; i < 3; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    Vector operator/(double scalar) const {
        Vector ans = *this;
        ans /= scalar;
        return ans;
    }

    bool operator==(const Vector& other) const {
        for (size_t i = 0; i < 3; ++i) {
            if (std::abs(data_[i] - other.data_[i]) >= kEps) {
                return false;
            }
        }
        return true;
    }

private:
    std::array<double, 3> data_;
};

inline double DotProduct(const Vector& lhs, const Vector& rhs) {
    double ans = 0;
    for (size_t i = 0; i < 3; ++i) {
        ans += lhs[i] * rhs[i];
    }
    return ans;
}

inline Vector CrossProduct(const Vector& a, const Vector& b) {
    Vector ans;
    ans[0] = a[1] * b[2] - a[2] * b[1];
    ans[1] = a[2] * b[0] - a[0] * b[2];
    ans[2] = a[0] * b[1] - a[1] * b[0];
    return ans;
}

inline double Length(const Vector& vec) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

inline double Dist(const Vector& a, const Vector& b) {
    return Length(a - b);
}