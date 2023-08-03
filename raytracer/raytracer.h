#pragma once

#include <image.h>
#include <camera_options.h>
#include <render_options.h>
#include <vector.h>
#include <ray.h>
#include <scene.h>
#include <intersection.h>
#include <geometry.h>

#include <string>
#include <cmath>
#include <array>
#include <optional>
#include <algorithm>

const double kAlign = 0.001;

std::array<std::array<double, 4>, 4> GenMatrix(const Vector& from, const Vector& to) {
    std::array<std::array<double, 4>, 4> m;
    Vector forward = from - to;
    forward.Normalize();
    Vector up({0, 1, 0});
    if (up == forward) {
        up = Vector({0, 0, -1});
    }
    if (forward == Vector({0, -1, 0})) {
        up = Vector({0, 0, 1});
    }
    Vector right = CrossProduct(up, forward);
    right.Normalize();
    Vector newup = CrossProduct(forward, right);

    m[0][0] = right[0], m[0][1] = right[1], m[0][2] = right[2], m[0][3] = 0;
    m[1][0] = newup[0], m[1][1] = newup[1], m[1][2] = newup[2], m[1][3] = 0;
    m[2][0] = forward[0], m[2][1] = forward[1], m[2][2] = forward[2], m[2][3] = 0;
    m[3][0] = from[0], m[3][1] = from[1], m[3][2] = from[2], m[3][3] = 1;
    return m;
}

Vector TransformToWorld(const std::array<std::array<double, 4>, 4>& m, const Vector& v) {
    double a, b, c, w;
    a = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0] + m[3][0];
    b = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1] + m[3][1];
    c = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2] + m[3][2];
    w = v[0] * m[0][3] + v[1] * m[1][3] + v[2] * m[2][3] + m[3][3];
    Vector result({a, b, c});
    result /= w;
    return result;
}

Vector TraceRay(const Ray& ray, const Scene& scene, int rem_depth, bool inside_object, const RenderMode& mode) {
    if (rem_depth == -1) {
        return {};
    }

    std::optional<Intersection> closest_intersection;
    std::optional<Object> closest_obj;
    std::optional<SphereObject> closest_sphere;
    double minimal_distance = 1000000000000000;

    const auto& objs = scene.GetObjects();
    for (const auto& elem : objs) {
        const auto& tmp = GetIntersection(ray, elem.polygon);
        if (!tmp.has_value()) {
            continue;
        }
        if (tmp->GetDistance() > kEps && tmp->GetDistance() < minimal_distance - kEps) {
            minimal_distance = tmp->GetDistance();
            closest_intersection = tmp;
            closest_obj = elem;
        }
    }

    const auto& sphere_objs = scene.GetSphereObjects();
    for (const auto& elem : sphere_objs) {
        const auto& tmp = GetIntersection(ray, elem.sphere);
        if (!tmp.has_value()) {
            continue;
        }
        if (tmp->GetDistance() > kEps && tmp->GetDistance() < minimal_distance - kEps) {
            minimal_distance = tmp->GetDistance();
            closest_sphere = elem;
            closest_intersection = tmp;
            closest_obj = std::nullopt;
        }
    }

    //  No intersection
    if (!closest_obj.has_value() && !closest_sphere.has_value()) {
        if (mode == RenderMode::kDepth) {
            return {-1, -1, -1};
        }
        if (mode == RenderMode::kNormal) {
            return {-228, -228, -228};
        }
        return {};
    }

    Vector result;
    const auto& lights = scene.GetLights();
    Intersection inter = closest_intersection.value();
    const Material* material = nullptr;

    if (closest_obj.has_value()) {
        //  Intersection with triangle
        material = closest_obj->material;

        //  Using barycentric coords
        if (closest_obj->GetNormal(0) != nullptr && closest_obj->GetNormal(1) != nullptr &&
            closest_obj->GetNormal(2) != nullptr) {
            Vector bars = GetBarycentricCoords(closest_obj->polygon, inter.GetPosition());
            inter.SetNormal(*closest_obj->GetNormal(0) * bars[0] +
                            *closest_obj->GetNormal(1) * bars[1] +
                            *closest_obj->GetNormal(2) * bars[2]);
        }
    } else {
        //  Intersection with sphere
        material = closest_sphere->material;
    }

    if (mode == RenderMode::kDepth) {
        return {inter.GetDistance(), inter.GetDistance(), inter.GetDistance()};
    }
    if (mode == RenderMode::kNormal) {
        return inter.GetNormal();
    }

    //  Ka + Ke
    result += material->ambient_color;
    result += material->intensity;

    Vector inter_aligned = inter.GetPosition() + inter.GetNormal() * kAlign * 2.0;
    Vector lights_sum;
    if (material->albedo[0] > kEps) {
        for (const auto& light : lights) {
            bool is_shadow = false;
            Ray light_ray(light.position, inter_aligned - light.position);
            double min_dist = Dist(inter.GetPosition(), light.position);

            for (const auto& obj : objs) {
                const auto& tmp = GetIntersection(light_ray, obj.polygon);
                if (!tmp.has_value()) {
                    continue;
                }
                if (tmp->GetDistance() < min_dist - kEps) {
                    is_shadow = true;
                    break;
                }
            }

            if (is_shadow) {
                continue;
            }

            for (const auto& sphere : sphere_objs) {
                const auto& tmp = GetIntersection(light_ray, sphere.sphere);
                if (!tmp.has_value()) {
                    continue;
                }
                if (tmp->GetDistance() < min_dist - kEps) {
                    is_shadow = true;
                    break;
                }
            }

            if (is_shadow) {
                continue;
            }

            Vector vl = light_ray.GetDirection() * (-1);
            vl.Normalize();
            //  Diffuse
            lights_sum += light.intensity * material->diffuse_color *
                          std::max(DotProduct(inter.GetNormal(), vl), static_cast<double>(0));

            Vector ve, vr;
            vr = inter.GetNormal() * DotProduct(inter.GetNormal(), vl) * 2.0 - vl;
            ve = ray.GetDirection() * (-1);
            ve.Normalize();
            //  Specular
            lights_sum += light.intensity * material->specular_color *
                          std::pow(std::max(DotProduct(ve, vr), static_cast<double>(0)),
                                   material->specular_exponent);
        }

        result += lights_sum * material->albedo[0];
    }

    if (!inside_object) {
        //  Reflect
        if (material->albedo[1] > kEps) {
            Ray reflected_ray(inter_aligned, Reflect(ray.GetDirection(), inter.GetNormal()));
            result +=
                TraceRay(reflected_ray, scene, rem_depth - 1, inside_object, mode) * material->albedo[1];
        }
    }

    //  Refract
    double eta = (inside_object ? material->refraction_index : 1.0 / material->refraction_index);
    double tr = (inside_object ? 1.0 : material->albedo[2]);
    if (tr > kEps) {
        std::optional<Vector> refrac_vec = Refract(ray.GetDirection(), inter.GetNormal(), eta);
        if (refrac_vec.has_value()) {
            Ray refracted_ray(inter_aligned - inter.GetNormal() * kAlign * 4.0, refrac_vec.value());
            result += TraceRay(refracted_ray, scene, rem_depth - 1,
                               (inside_object ^ closest_sphere.has_value()), mode) *
                      tr;
        }
    }
    return result;
}

Image Render(const std::string& filename, const CameraOptions& camera_options,
             const RenderOptions& render_options) {
    Image result(camera_options.screen_width, camera_options.screen_height);
    std::vector<std::vector<Vector>> screen(camera_options.screen_height,
                                            std::vector<Vector>(camera_options.screen_width));

    auto camera_to_world = GenMatrix(camera_options.look_from, camera_options.look_to);
    const auto& scene = ReadScene(filename);
    /*const auto& ppp = scene.GetSphereObjects();
    for (const auto& elem : ppp) {
        for (size_t i = 0; i < 3; ++i) {
            std::cout << elem.material->albedo[i] << " ";
        }
        std::cout << "\n";
    }*/

    double scale = std::tan(camera_options.fov / 2.0);
    double image_aspect_ratio =
        static_cast<double>(camera_options.screen_width) / camera_options.screen_height;
    Vector camera_origin = TransformToWorld(camera_to_world, Vector({0, 0, 0}));
    double max_value = -1000000000000000;

    /*for (int pixel_x = 400; pixel_x <= 400; ++pixel_x) {
        for (int pixel_y = 300; pixel_y <= 300; ++pixel_y) {*/
    for (int pixel_x = 0; pixel_x < camera_options.screen_width; ++pixel_x) {
        for (int pixel_y = 0; pixel_y < camera_options.screen_height; ++pixel_y) {
            double pixel_ndc_x = (static_cast<double>(pixel_x) + 0.5) / camera_options.screen_width;
            double pixel_ndc_y =
                (static_cast<double>(pixel_y) + 0.5) / camera_options.screen_height;
            double pixel_screen_x = pixel_ndc_x * 2 - 1.0;
            double pixel_screen_y = 1.0 - pixel_ndc_y * 2;
            double pixel_camera_x = pixel_screen_x * image_aspect_ratio * scale;
            double pixel_camera_y = pixel_screen_y * scale;
            Vector pixel =
                TransformToWorld(camera_to_world, Vector({pixel_camera_x, pixel_camera_y, -1}));
            //  trace ray and fill image
            screen[pixel_y][pixel_x] = TraceRay(Ray(camera_origin, pixel - camera_origin), scene,
                                                render_options.depth - 1, false, render_options.mode);
            max_value = std::max({max_value, screen[pixel_y][pixel_x][0],
                                  screen[pixel_y][pixel_x][1], screen[pixel_y][pixel_x][2]});
        }
    }

    //  Postrprocessing
    for (int pixel_x = 0; pixel_x < camera_options.screen_width; ++pixel_x) {
        for (int pixel_y = 0; pixel_y < camera_options.screen_height; ++pixel_y) {
            Vector& vec = screen[pixel_y][pixel_x];
            if (render_options.mode == RenderMode::kDepth && vec == Vector({-1, -1, -1})) {
                vec = Vector({1, 1, 1});
            } else {
                if (render_options.mode == RenderMode::kNormal && vec == Vector({-228, -228, -228})) {
                    vec = Vector();
                } else {
                    for (size_t i = 0; i < 3; ++i) {
                        if (render_options.mode == RenderMode::kDepth) {
                            vec[i] /= max_value;
                        } else {
                            if (render_options.mode == RenderMode::kNormal) {
                                vec[i] = (0.5 * vec[i] + 0.5);
                            } else {
                                vec[i] = (vec[i] * (1.0 + (vec[i] / (max_value * max_value)))) / (vec[i] + 1.0);
                                vec[i] = std::pow(vec[i], 1.0 / 2.2);
                            }
                        }
                    }
                }
            }
            RGB tmp;
            tmp.r = static_cast<int>(vec[0] * 255);
            tmp.g = static_cast<int>(vec[1] * 255);
            tmp.b = static_cast<int>(vec[2] * 255);
            result.SetPixel(tmp, pixel_y, pixel_x);
        }
    }

    return result;
}
