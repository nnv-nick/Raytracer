#pragma once

#include <vector.h>

class Triangle {
public:
    Triangle(const std::initializer_list<Vector>& list) {
        size_t ind = 0;
        for (const auto& elem : list) {
            vertices_[ind++] = elem;
        }
    }
    double Area() const {
        double a = Length(vertices_[0] - vertices_[1]);
        double b = Length(vertices_[1] - vertices_[2]);
        double c = Length(vertices_[2] - vertices_[0]);
        double p = (a + b + c) / 2.0;
        return sqrt(p * (p - a) * (p - b) * (p - c));
    }

    const Vector& GetVertex(size_t ind) const {
        assert(ind < 3);
        return vertices_[ind];
    }

private:
    std::array<Vector, 3> vertices_;
};
