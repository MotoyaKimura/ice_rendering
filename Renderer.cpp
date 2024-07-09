//
// Created by sty-g on 2023/11/16.
//

#include "Renderer.h"
#include "TriMesh.h"
#include <iostream>

Renderer::Renderer(const RTCScene &scene, const Object &objects, Camera camera, Color bgColor)
        : scene(scene), objects(objects), camera(std::move(camera)), bgColor(std::move(bgColor)), engine(0), dist(0, 1) {
}

double Renderer::rand() const {
    return dist(engine);
}

Image Renderer::render() const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    std::cout << "width: " << image.width << ", height: " << image.height << std::endl;

    RTCIntersectArguments arguments;
    rtcInitIntersectArguments(&arguments);

    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Color color;
            RTCRayHit rayHit;
            camera.filmView(p_x, p_y, rayHit);
            rayHit.ray.tnear = 1e-3f;
            rayHit.ray.time = 0.0f;
            rayHit.ray.tfar = std::numeric_limits<float>::max();
            rayHit.ray.mask = 0xFFFFFFFF;
            rayHit.ray.flags = 0;
            rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rayHit.hit.primID = RTC_INVALID_GEOMETRY_ID;
            rtcIntersect1(scene, &rayHit, &arguments);

            if(rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                int sum = 0;
                for(int i = 0; i < objects.meshes.size(); i++) {
                    if(rayHit.hit.primID < objects.meshes[i].triangles.size() + sum) {
                        color.x() = objects.meshes.at(i).material.Kd.x();
                        color.y() = objects.meshes.at(i).material.Kd.y();
                        color.z() = objects.meshes.at(i).material.Kd.z();
                        break;
                    } else {
                        sum += objects.meshes[i].triangles.size();
                    }
                }
            } else {
                color = bgColor;
            }
            image.pixels[p_idx] = color;
        }
    }

    return image;
}

Image Renderer::directIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    std::cout << "width: " << image.width << ", height: " << image.height << std::endl;

#pragma omp parallel for
    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            for(int i = 0; i < samples; ++i) {
//                std::cout << "x: " << p_x << ", y: " << p_y << ", i: " << i << std::endl;
                RTCRayHit rayHit;
                camera.filmView(p_x, p_y, rayHit);
                Color reflectRadiance = Color::Zero();
                reflectRadiance = russian(rayHit);
                image.pixels[p_idx] += reflectRadiance;
            }
            image.pixels[p_idx] /= static_cast<double>(samples);
        }
    }
    return image;
}

Image Renderer::normalRender() const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());

    RTCIntersectArguments arguments;
    rtcInitIntersectArguments(&arguments);

    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Color color;
            RTCRayHit rayHit;
            camera.filmView(p_x, p_y, rayHit);
            rayHit.ray.tnear = 1e-3f;
            rayHit.ray.time = 0.0f;
            rayHit.ray.tfar = std::numeric_limits<float>::max();
            rayHit.ray.mask = 0xFFFFFFFF;
            rayHit.ray.flags = 0;
            rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rayHit.hit.primID = RTC_INVALID_GEOMETRY_ID;
            rtcIntersect1(scene, &rayHit, &arguments);

            Eigen::Vector3d in_ray_dir = Eigen::Vector3d(rayHit.ray.dir_x, rayHit.ray.dir_y, rayHit.ray.dir_z).normalized();
            Eigen::Vector3d in_ray_normal = Eigen::Vector3d(rayHit.hit.Ng_x, rayHit.hit.Ng_y, rayHit.hit.Ng_z).normalized();

            if(in_ray_normal.dot(in_ray_dir) > 0) {
                in_ray_normal = -in_ray_normal;
            }

            if(rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                color.x() = in_ray_normal.x() / 2 + 0.5;
                color.y() = in_ray_normal.y() / 2 + 0.5;
                color.z() = in_ray_normal.z() / 2 + 0.5;
            } else {
                color = bgColor;
            }
            image.pixels[p_idx] = color;
        }
    }

    return image;
}

//void Renderer::diffuseSample(RTCRayHit in_ray, RTCRayHit &out_ray, const TriMesh &triMesh, const int triangle_num, const int vert_sum) const {
void Renderer::diffuseSample(RTCRayHit in_ray, RTCRayHit &out_ray) const {
    const double phi = 2.0 * EIGEN_PI * rand();
    const double theta = asin(sqrt(rand()));

    const double _x = sin(theta) * cos(phi);
    const double _y = cos(theta);
    const double _z = sin(theta) * sin(phi);

    Eigen::Vector3d u, v;
    const Eigen::Vector3d in_ray_dir = Eigen::Vector3d(in_ray.ray.dir_x, in_ray.ray.dir_y, in_ray.ray.dir_z).normalized();
    Eigen::Vector3d in_ray_normal = Eigen::Vector3d(in_ray.hit.Ng_x, in_ray.hit.Ng_y, in_ray.hit.Ng_z).normalized();

    if(in_ray_normal.dot(in_ray_dir) > 0) {
        in_ray_normal = -in_ray_normal;
    }

    computeLocalFrame(in_ray_normal, u, v);
    const Eigen::Vector3d dir = _x * u + _y * in_ray_normal + _z * v;
    out_ray.ray.dir_x = dir.x();
    out_ray.ray.dir_y = dir.y();
    out_ray.ray.dir_z = dir.z();

    out_ray.ray.org_x = in_ray.ray.org_x + in_ray.ray.dir_x * in_ray.ray.tfar + dir.x() * 1e-3;
    out_ray.ray.org_y = in_ray.ray.org_y + in_ray.ray.dir_y * in_ray.ray.tfar + dir.y() * 1e-3;
    out_ray.ray.org_z = in_ray.ray.org_z + in_ray.ray.dir_z * in_ray.ray.tfar + dir.z() * 1e-3;
}

void Renderer::computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v) {
    if(fabs(w.x()) > 1e-3)
        u = Eigen::Vector3d ::UnitY().cross(w).normalized();
    else
        u = Eigen::Vector3d::UnitX().cross(w).normalized();

    v = w.cross(u);
}

void Renderer::specularReflection(const RTCRayHit &in_ray, RTCRayHit &out_ray) const {
    const Eigen::Vector3d in_ray_dir = Eigen::Vector3d(in_ray.ray.dir_x, in_ray.ray.dir_y, in_ray.ray.dir_z).normalized();
    Eigen::Vector3d in_ray_normal = Eigen::Vector3d(in_ray.hit.Ng_x, in_ray.hit.Ng_y, in_ray.hit.Ng_z).normalized();

    if(in_ray_normal.dot(in_ray_dir) > 0) {
        in_ray_normal = -in_ray_normal;
    }

    const Eigen::Vector3d dir = in_ray_dir - 2 * in_ray_dir.dot(in_ray_normal) * in_ray_normal;
    out_ray.ray.dir_x = dir.x();
    out_ray.ray.dir_y = dir.y();
    out_ray.ray.dir_z = dir.z();

    out_ray.ray.org_x = in_ray.ray.org_x + in_ray.ray.dir_x * in_ray.ray.tfar + dir.x() * 1e-3;
    out_ray.ray.org_y = in_ray.ray.org_y + in_ray.ray.dir_y * in_ray.ray.tfar + dir.y() * 1e-3;
    out_ray.ray.org_z = in_ray.ray.org_z + in_ray.ray.dir_z * in_ray.ray.tfar + dir.z() * 1e-3;
}

void Renderer::calc_refract(const RTCRayHit &in_ray, RTCRayHit &out_ray) const {
    const Eigen::Vector3d in_ray_dir = Eigen::Vector3d(in_ray.ray.dir_x, in_ray.ray.dir_y, in_ray.ray.dir_z).normalized();
    const Eigen::Vector3d in_ray_normal = Eigen::Vector3d(in_ray.hit.Ng_x, in_ray.hit.Ng_y, in_ray.hit.Ng_z).normalized();

    const double ref_idx = 1.3;
    double ni_over_nt;
    double reflect_prob;
    double cosine = -in_ray_dir.dot(in_ray_normal);

    if(cosine < 0) {
        cosine = -in_ray_dir.dot(-in_ray_normal);
        ni_over_nt = ref_idx;
    } else {
        ni_over_nt = 1 / ref_idx;
    }

    if(refract(in_ray_dir, in_ray_normal, ni_over_nt, out_ray)) {
        reflect_prob = schlick(cosine, ref_idx);
    } else {
        reflect_prob = 1.0;
    }

    if((double)rand() < reflect_prob) {
        specularReflection(in_ray, out_ray);
    } else {
        out_ray.ray.org_x = in_ray.ray.org_x + in_ray.ray.dir_x * in_ray.ray.tfar + out_ray.ray.dir_x * 1e-3;
        out_ray.ray.org_y = in_ray.ray.org_y + in_ray.ray.dir_y * in_ray.ray.tfar + out_ray.ray.dir_y * 1e-3;
        out_ray.ray.org_z = in_ray.ray.org_z + in_ray.ray.dir_z * in_ray.ray.tfar + out_ray.ray.dir_z * 1e-3;
    }
}

bool Renderer::refract(const Eigen::Vector3d &dir, Eigen::Vector3d normal, const double ni_over_nt, RTCRayHit &out_ray) const {

    if(-dir.dot(normal) < 0) {
        normal = -normal;
    }

    const double dt = -dir.dot(normal);
    const double discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);

    if(discriminant > 0) {
        const Eigen::Vector3d out_ray_dir = ni_over_nt * (dir + dt * normal) - sqrt(discriminant) * normal;
        out_ray.ray.dir_x = out_ray_dir.x();
        out_ray.ray.dir_y = out_ray_dir.y();
        out_ray.ray.dir_z = out_ray_dir.z();
        return true;
    } else {
        return false;
    }
}

double Renderer::schlick(const double cosine, const double ref_idx) const {
    const double r0 = pow((1 - ref_idx) / (1 + ref_idx), 2);
    return r0 + (1 + r0) * pow((1 - cosine), 5);
}

Color Renderer::russian(RTCRayHit rayHit) const {
    RTCRayHit _rayHit;
    const double rand01 = rand();
    int obj_idx = 0;

    rayHit.ray.tnear = 1e-3f;
    rayHit.ray.time = 0.0f;
    rayHit.ray.tfar = std::numeric_limits<float>::max();
    rayHit.ray.mask = 0xFFFFFFFF;
    rayHit.ray.flags = 0;
    rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayHit.hit.primID = RTC_INVALID_GEOMETRY_ID;

    RTCIntersectArguments arguments;
    rtcInitIntersectArguments(&arguments);
    rtcIntersect1(scene, &rayHit, &arguments);

    _rayHit.ray.tnear = 1e-3f;
    _rayHit.ray.time = 0.0f;
    _rayHit.ray.tfar = std::numeric_limits<float>::max();
    _rayHit.ray.mask = 0xFFFFFFFF;
    _rayHit.ray.flags = 0;
    _rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    _rayHit.hit.primID = RTC_INVALID_GEOMETRY_ID;

    if(rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {

        int tri_sum = 0;
        for(int i = 0; i < objects.meshes.size(); i++) {
            if(rayHit.hit.primID < objects.meshes.at(i).triangles.size() + tri_sum) {
                obj_idx = i;
                break;
            } else {
                tri_sum += objects.meshes.at(i).triangles.size();
            }
        }

        const double kd_max = objects.meshes.at(obj_idx).material.Kd.maxCoeff();
        const double ks_max = objects.meshes.at(obj_idx).material.Ks.maxCoeff();
        const double kt = objects.meshes.at(obj_idx).material.d;    /// アルファの値　デフォルト；1　だから、1-[アルファ]に直す

        Color reflection = Color::Zero();
        double probability;

        if(rand01 < kd_max) {
            probability = kd_max;
            if(objects.meshes.at(obj_idx).material.name != "light") {
                diffuseSample(rayHit, _rayHit);
            } else {
                _rayHit.ray.dir_x = rayHit.ray.dir_x;
                _rayHit.ray.dir_y = rayHit.ray.dir_y;
                _rayHit.ray.dir_z = rayHit.ray.dir_z;
                _rayHit.ray.org_x = rayHit.ray.org_x + rayHit.ray.dir_x * rayHit.ray.tfar + rayHit.ray.dir_x * 1e-3;
                _rayHit.ray.org_y = rayHit.ray.org_y + rayHit.ray.dir_y * rayHit.ray.tfar + rayHit.ray.dir_y * 1e-3;
                _rayHit.ray.org_z = rayHit.ray.org_z + rayHit.ray.dir_z * rayHit.ray.tfar + rayHit.ray.dir_z * 1e-3;
            }
//            diffuseSample(rayHit, _rayHit);
            const auto incomingLight = russian(_rayHit);
            reflection = objects.meshes.at(obj_idx).material.Kd.cwiseProduct(incomingLight) / probability;
        } else if(rand01 < kd_max + ks_max) {
            probability = ks_max;
            specularReflection(rayHit, _rayHit);
            const auto incomingLight = russian(_rayHit);
            reflection = objects.meshes.at(obj_idx).material.Ks.cwiseProduct(incomingLight) / probability;
        } else if(rand01 < kd_max + ks_max + kt) {
//            std::cout << objects.meshes.at(obj_idx).material.name << std::endl;
            probability = kt;
            calc_refract(rayHit, _rayHit);
            const auto incomingLight = russian(_rayHit);
            reflection = objects.meshes.at(obj_idx).material.d * incomingLight / probability;
        }
        const Color result = objects.meshes.at(obj_idx).material.Ke + reflection;
        return result;
    }
    return bgColor;
}