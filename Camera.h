//
// Created by sty-g on 2023/11/16.
//

#ifndef HELLO_EMBREE_CAMERA_H
#define HELLO_EMBREE_CAMERA_H

#include "Film.h"
#include <embree4/rtcore.h>
#include <embree4/rtcore_ray.h>
#include <embree4/rtcore_scene.h>

class Camera {
private:
    double focalLength;
    Eigen::Vector3d org, right, up;
    Film film;

public:
    Camera(Eigen::Vector3d position, const Eigen::Vector3d &dir, const int &resolutionHeight, double aspectRatio, double verticalFov, double focalLength=1.0);
    void filmView(const unsigned int &p_x, const unsigned int &p_y, RTCRayHit &rayHit) const;

    const Film &getFilm() const;

    Eigen::Vector3d dir;
};


#endif //HELLO_EMBREE_CAMERA_H
