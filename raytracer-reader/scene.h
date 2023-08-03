#pragma once

#include <material.h>
#include <vector.h>
#include <object.h>
#include <light.h>

#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <cassert>
#include <utility>
#include <string>

class Scene {
public:
    Scene() {
    }

    const std::vector<Object>& GetObjects() const {
        return objs_;
    }

    const std::vector<SphereObject>& GetSphereObjects() const {
        return spheres_;
    }

    const std::vector<Light>& GetLights() const {
        return lights_;
    }

    const std::map<std::string, Material>& GetMaterials() const {
        return name_to_mat_;
    }

    void SetMaterial(const std::string& name, const Material& mat) {
        name_to_mat_[name] = mat;
    }

    const Material* GetMaterial(const std::string& name) const {
        assert(name_to_mat_.count(name));
        return &name_to_mat_.at(name);
    }

    void SetObjs(const std::vector<Object>& objs) {
        objs_ = objs;
    }

    void SetSpheres(const std::vector<SphereObject>& spheres) {
        spheres_ = spheres;
    }

    void SetLights(const std::vector<Light>& lights) {
        lights_ = lights;
    }

private:
    std::vector<Object> objs_;
    std::vector<SphereObject> spheres_;
    std::vector<Light> lights_;
    std::map<std::string, Material> name_to_mat_;
};

inline std::map<std::string, Material> ReadMaterials(std::string_view filename) {
    std::ifstream fin(filename.data());
    assert(fin.is_open());
    std::map<std::string, Material> result;
    std::string s, cur_name;
    Material mat, emp;
    double x, y, z;
    while (!fin.eof()) {
        fin >> s;
        if (s == "newmtl") {
            if (!cur_name.empty()) {
                result[cur_name] = mat;
            }
            fin >> cur_name;
            mat = emp;
            mat.name = cur_name;
            continue;
        }
        if (s == "Ka") {
            fin >> x >> y >> z;
            mat.ambient_color = Vector({x, y, z});
            continue;
        }
        if (s == "Kd") {
            fin >> x >> y >> z;
            mat.diffuse_color = Vector({x, y, z});
            continue;
        }
        if (s == "Ks") {
            fin >> x >> y >> z;
            mat.specular_color = Vector({x, y, z});
            continue;
        }
        if (s == "Ke") {
            fin >> x >> y >> z;
            mat.intensity = Vector({x, y, z});
            continue;
        }
        if (s == "Ns") {
            fin >> x;
            mat.specular_exponent = x;
            continue;
        }
        if (s == "Ni") {
            fin >> x;
            mat.refraction_index = x;
            continue;
        }
        if (s == "al") {
            fin >> x >> y >> z;
            mat.albedo = std::array<double, 3>({x, y, z});
            continue;
        }
        getline(fin, s);
    }
    if (!cur_name.empty()) {
        result[cur_name] = mat;
    }
    return result;
}

inline std::vector<std::pair<int, int>> ParseFace(std::string& s) {
    std::vector<std::pair<int, int>> result;
    std::string word;
    s += ' ';
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '/' || s[i] == '-' || (s[i] >= '0' && s[i] <= '9')) {
            word += s[i];
        } else {
            if (!word.empty()) {
                int pos1 = -1, pos2 = -1;
                for (int j = 0; j < static_cast<int>(word.size()); ++j) {
                    if (word[j] != '/') {
                        continue;
                    }
                    if (pos1 == -1) {
                        pos1 = j;
                    }
                    pos2 = j;
                }
                if (pos1 == -1) {
                    result.push_back(std::make_pair(std::stoi(word), 0));
                } else {
                    std::pair<int, int> tmp;
                    tmp.first = std::stoi(word.substr(0, pos1));
                    int sz = static_cast<int>(word.size());
                    if (pos2 == sz - 1) {
                        tmp.second = 0;
                    } else {
                        tmp.second = std::stoi(word.substr(pos2 + 1, sz - pos2 - 1));
                    }
                    result.push_back(tmp);
                }
            }
            word = "";
        }
    }
    return result;
}

inline int ChInd(int x, int p) {
    if (x > 0) {
        return x - 1;
    }
    return p + x;
}

inline Scene ReadScene(const std::string& filename) {
    std::string file_pref = filename.substr(0, filename.rfind('/') + 1);
    std::ifstream fin(filename.data());
    assert(fin.is_open());
    Scene result;
    std::string s;
    std::vector<Object> objs;
    std::vector<SphereObject> spheres;
    std::vector<Light> lights;
    std::vector<Vector> vertices, normals;
    const Material* cur_mat = nullptr;
    double x, y, z;
    while (!fin.eof()) {
        fin >> s;
        if (s == "v") {
            fin >> x >> y >> z;
            vertices.push_back({x, y, z});
            continue;
        }
        if (s == "vn") {
            fin >> x >> y >> z;
            Vector tmp({x, y, z});
            tmp.Normalize();
            normals.push_back(tmp);
            continue;
        }
        if (s == "f") {
            getline(fin, s);
            std::vector<std::pair<int, int>> tmp = ParseFace(s);
            for (size_t i = 1; i < tmp.size() - 1; ++i) {
                int sz = static_cast<int>(vertices.size());
                objs.emplace_back(
                    cur_mat,
                    Triangle({vertices[ChInd(tmp[0].first, sz)], vertices[ChInd(tmp[i].first, sz)],
                              vertices[ChInd(tmp[i + 1].first, sz)]}),
                    std::unordered_map<size_t, Vector>({}));
                sz = static_cast<int>(normals.size());
                if (tmp[0].second != 0) {
                    objs.back().ind_to_norm[0] = normals[ChInd(tmp[0].second, sz)];
                }
                if (tmp[i].second != 0) {
                    objs.back().ind_to_norm[1] = normals[ChInd(tmp[i].second, sz)];
                }
                if (tmp[i + 1].second != 0) {
                    objs.back().ind_to_norm[2] = normals[ChInd(tmp[i + 1].second, sz)];
                }
            }
            continue;
        }
        if (s == "mtllib") {
            fin >> s;
            auto tmp = ReadMaterials(file_pref + s);
            for (const auto& [key, value] : tmp) {
                result.SetMaterial(key, value);
            }
            continue;
        }
        if (s == "usemtl") {
            fin >> s;
            cur_mat = result.GetMaterial(s);
            continue;
        }
        if (s == "S") {
            double r;
            fin >> x >> y >> z >> r;
            spheres.emplace_back(cur_mat, Sphere(Vector({x, y, z}), r));
            continue;
        }
        if (s == "P") {
            double r, g, b;
            fin >> x >> y >> z >> r >> g >> b;
            lights.emplace_back(Vector({x, y, z}), Vector({r, g, b}));
            continue;
        }
        getline(fin, s);
    }
    result.SetLights(lights);
    result.SetObjs(objs);
    result.SetSpheres(spheres);
    return result;
}
