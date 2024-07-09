//
// Created by sty-g on 2023/11/20.
//

#ifndef HELLO_EMBREE_TRIMESH_H
#define HELLO_EMBREE_TRIMESH_H

#include <Eigen/Dense>

struct Material {
    double Ns;
    Eigen::Vector3d Ka;
    Eigen::Vector3d Kd;
    Eigen::Vector3d Ks;
    Eigen::Vector3d Ke;
    double Ni;
    double d;
    int illum;
    std::string name;
};

struct TriMesh {
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertices;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertex_normals;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> tex_coords;
    std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> triangles;
//    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertex_colors;
    Material material;
};

struct Object {
    std::vector<TriMesh> meshes;
};

#endif //HELLO_EMBREE_TRIMESH_H
