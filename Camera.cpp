//
// Created by sty-g on 2023/11/16.
//

#include "Camera.h"
#include <iostream>

Camera::Camera(Eigen::Vector3d position, const Eigen::Vector3d &dir, const int &resolutionHeight, double aspectRatio, double verticalFov, double focalLength)
    : org(std::move(position)), dir(dir.normalized()), focalLength(focalLength) {
    const auto theta = verticalFov * EIGEN_PI / 180.0;

    const double filmHeight = 2.0 * tan(theta / 2.0) * focalLength;
    film = Film(resolutionHeight, aspectRatio, filmHeight);

    right = this->dir.cross(Eigen::Vector3d::UnitY());
    up = right.cross(this->dir).normalized();
}

void Camera::filmView(const unsigned int &p_x, const unsigned int &p_y, RTCRayHit &rayHit) const {
    const auto pixelLocalPos = film.pixelLocalPosition(p_x, p_y);
    const Eigen::Vector3d pixelWorldPos = org + film.filmSize.x() * right * (pixelLocalPos.x() - 0.5) + film.filmSize.y() * up * (0.5 - pixelLocalPos.y()) + focalLength * dir;
    rayHit.ray.org_x = pixelWorldPos.x();
    rayHit.ray.org_y = pixelWorldPos.y();
    rayHit.ray.org_z = pixelWorldPos.z();
    rayHit.ray.dir_x = (pixelWorldPos.x() - org.x()) / sqrt(pow(pixelWorldPos.x() - org.x(), 2) + pow(pixelWorldPos.y() - org.y(), 2) + pow(pixelWorldPos.z() - org.z(), 2));
    rayHit.ray.dir_y = (pixelWorldPos.y() - org.y()) / sqrt(pow(pixelWorldPos.x() - org.x(), 2) + pow(pixelWorldPos.y() - org.y(), 2) + pow(pixelWorldPos.z() - org.z(), 2));
    rayHit.ray.dir_z = (pixelWorldPos.z() - org.z()) / sqrt(pow(pixelWorldPos.x() - org.x(), 2) + pow(pixelWorldPos.y() - org.y(), 2) + pow(pixelWorldPos.z() - org.z(), 2));
}

const Film &Camera::getFilm() const {
    return film;
}