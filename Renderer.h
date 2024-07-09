//
// Created by sty-g on 2023/11/16.
//

#ifndef HELLO_EMBREE_RENDERER_H
#define HELLO_EMBREE_RENDERER_H

#include <vector>
#include "Camera.h"
#include "TriMesh.h"
#include <random>

class Renderer {
public:
    RTCScene scene;
    Object objects;

    Camera camera;
    Color bgColor;

    mutable std::mt19937_64 engine;
    mutable std::uniform_real_distribution<> dist;

    Renderer(const RTCScene &scene, const Object &objects, Camera camera, Color bgColor=Color::Zero());

    double rand() const;

    Image render() const;

    Image directIlluminationRender(const unsigned int &samples) const;

    Image normalRender() const;

    void diffuseSample(RTCRayHit in_ray, RTCRayHit &out_ray) const;

    static void computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v);

    void specularReflection(const RTCRayHit &in_ray, RTCRayHit &out_ray) const;

    void calc_refract(const RTCRayHit &in_ray, RTCRayHit &out_ray) const;

    bool refract(const Eigen::Vector3d &dir, Eigen::Vector3d normal, const double ni_over_nt, RTCRayHit &out_ray) const;

    double schlick(const double cosine, const double ref_idx) const;

    Color russian(RTCRayHit rayHit) const;
};


#endif //HELLO_EMBREE_RENDERER_H
