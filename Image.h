//
// Created by sty-g on 2023/11/16.
//

#ifndef HELLO_EMBREE_IMAGE_H
#define HELLO_EMBREE_IMAGE_H

#include "Color.h"
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

class Image {
public:
    int width;
    int height;

    Color *pixels;

    Image() = default;

    Image(const int &width, const int &height);

    ~Image();

    void save(const std::string &fname) const;

    Image apply_gamma_correction() const;

    cv::Mat toCvMat() const;

    cv::Mat toCvMat8UC3() const;
};


#endif //HELLO_EMBREE_IMAGE_H
