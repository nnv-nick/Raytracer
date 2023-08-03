#pragma once

#include <triangle.h>
#include <material.h>
#include <sphere.h>

#include <unordered_map>

struct Object {
    const Material* material = nullptr;
    Triangle polygon;
    std::unordered_map<size_t, Vector> ind_to_norm;

    Object(const Material* ptr, const Triangle& tr, const std::unordered_map<size_t, Vector>& m1)
        : material(ptr), polygon(tr), ind_to_norm(m1) {
    }

    const Vector* GetNormal(size_t index) const {
        if (ind_to_norm.count(index)) {
            return &ind_to_norm.at(index);
        }
        return nullptr;
    }
};

struct SphereObject {
    const Material* material = nullptr;
    Sphere sphere;

    SphereObject(const Material* ptr, const Sphere& sp) : material(ptr), sphere(sp) {
    }
};
