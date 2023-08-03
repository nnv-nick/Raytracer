#pragma once

#include <vector.h>

class Intersection {
public:
    Intersection(Vector pos, Vector norm, double dist)
        : position_(pos), normal_(norm), distance_(dist) {
        normal_.Normalize();
    }

    const Vector& GetPosition() const {
        return position_;
    }

    const Vector& GetNormal() const {
        return normal_;
    }

    double GetDistance() const {
        return distance_;
    }

    void SetNormal(const Vector& norm){
        normal_ = norm;
        normal_.Normalize();
    }

private:
    Vector position_;
    Vector normal_;
    double distance_;
};
