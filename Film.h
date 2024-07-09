//
// Created by sty-g on 2023/11/16.
//

#ifndef HELLO_EMBREE_FILM_H
#define HELLO_EMBREE_FILM_H

#include <Eigen/Dense>
#include "Image.h"

struct Film {
    Eigen::Vector2d filmSize;
    Eigen::Vector2i resolution;

    Film() = default;

    Film(const int &resolutionHeight, const double &aspectRatio, const double &filmHeight)
    : resolution(Eigen::Vector2i{resolutionHeight * aspectRatio, resolutionHeight}),
    filmSize(Eigen::Vector2d{filmHeight * aspectRatio, filmHeight}) {}

    Eigen::Vector2d pixelLocalPosition(const unsigned int &x, const unsigned int &y) const {
        return Eigen::Vector2d{(x + 0.5) / resolution.x(), (y + 0.5) / resolution.y()};
    }
};

#endif //HELLO_EMBREE_FILM_H
