#include <vector.h>
#include <sphere.h>
#include <ray.h>
#include <intersection.h>
#include <triangle.h>

#include <optional>

std::optional<Intersection> GetIntersection(const Ray& ray, const Sphere& sphere) {
    Vector l = sphere.GetCenter() - ray.GetOrigin();

    //  Proverka na tupoy ugol
    Vector tmp = l;
    tmp.Normalize();
    if (DotProduct(tmp, ray.GetDirection()) < -kEps) {
        return std::nullopt;
    }

    double tc = DotProduct(l, ray.GetDirection());
    double d = sqrt(Length(l) * Length(l) - tc * tc);
    if (d - sphere.GetRadius() > kEps) {
        return std::nullopt;
    }
    double t1c = sqrt(sphere.GetRadius() * sphere.GetRadius() - d * d);
    double t1 = tc - t1c;
    if (t1 < -kEps) {
        t1 = tc + t1c;
    }
    Vector p = ray.GetOrigin() + ray.GetDirection() * t1;
    Vector norm = p - sphere.GetCenter();
    norm.Normalize();
    if (Dist(ray.GetOrigin(), sphere.GetCenter()) < sphere.GetRadius() - kEps) {
        norm *= -1;
    }
    return Intersection(p, norm, Dist(p, ray.GetOrigin()));
}

std::optional<Intersection> GetIntersection(const Ray& ray, const Triangle& triangle) {
    Vector vertex0 = triangle.GetVertex(0);
    Vector vertex1 = triangle.GetVertex(1);
    Vector vertex2 = triangle.GetVertex(2);
    Vector edge1, edge2, h, s, q;
    double a, f, u, v;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    h = CrossProduct(ray.GetDirection(), edge2);
    a = DotProduct(edge1, h);
    if (std::abs(a) < kEps) {
        return std::nullopt;
    }
    f = 1.0 / a;
    s = ray.GetOrigin() - vertex0;
    u = f * DotProduct(s, h);
    if (u < -kEps || u > 1.0 + kEps) {
        return std::nullopt;
    }
    q = CrossProduct(s, edge1);
    v = f * DotProduct(q, ray.GetDirection());
    if (v < -kEps || u + v > 1.0 + kEps) {
        return std::nullopt;
    }
    float t = f * DotProduct(q, edge2);
    if (t > kEps) {
        Vector p = ray.GetOrigin() + ray.GetDirection() * t;
        Vector norm = CrossProduct(edge1, edge2);
        norm.Normalize();
        if (DotProduct(ray.GetDirection(), norm) > 0) {
            norm *= -1;
        }
        return Intersection(p, norm, Dist(p, ray.GetOrigin()));
    }
    return std::nullopt;
}

std::optional<Vector> Refract(const Vector& ray, const Vector& normal, double eta) {
    double co = -DotProduct(ray, normal);
    double si = sqrt(1.0 - co * co);
    if (si * eta > 1.0 - kEps) {
        return std::nullopt;
    }
    return ray * eta + normal * (eta * co - sqrt(1.0 - eta * eta * (1.0 - co * co)));
}

Vector Reflect(const Vector& ray, const Vector& normal) {
    double co = -DotProduct(ray, normal);
    return ray + normal * co * 2;
}

Vector GetBarycentricCoords(const Triangle& triangle, const Vector& point) {
    Vector vertex0 = triangle.GetVertex(0);
    Vector vertex1 = triangle.GetVertex(1);
    Vector vertex2 = triangle.GetVertex(2);
    Triangle t1{vertex1, vertex2, point};
    Triangle t2{vertex2, vertex0, point};
    Triangle t3{vertex0, vertex1, point};
    Vector ans{t1.Area(), t2.Area(), t3.Area()};
    double tmp = ans[0] + ans[1] + ans[2];
    ans /= tmp;
    return ans;
}
