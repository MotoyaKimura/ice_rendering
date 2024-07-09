//
// Created by sty-g on 2023/11/16.
//

#include "Image.h"
#include <iostream>

Image::Image(const int &width, const int &height) : width(width), height(height) {
    pixels = new Color[width * height];
    for(int i = 0; i < width * height; i++) {
        pixels[i][0] = 0;
        pixels[i][1] = 0;
        pixels[i][2] = 0;
    }
}

Image::~Image() {
    delete[] pixels;
}

void Image::save(const std::string &fname) const {
    const auto mat = toCvMat8UC3();

    bool success = cv::imwrite(fname, mat);

    if(success) std::cout << "save a image is success." << std::endl;
    else std::cout << "save a image is failed." << std::endl;
}

Image Image::apply_gamma_correction() const {
    Image img(width, height);

    for(int i = 0; i < width * height; i++) {
        for(int j = 0; j < 3; j++) {
            if(pixels[i][j] <= 0.0031308) img.pixels[i][j] = pixels[i][j] * 12.92;
            else img.pixels[i][j] = pow(pixels[i][j], 1.0 / 2.4) * 1.055 - 0.055;
        }
    }

    return img;
}

cv::Mat Image::toCvMat() const {
    cv::Mat mat(height, width, CV_64FC3);

    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            auto& pixel = mat.at<cv::Vec3d>(i, j);
            Eigen::Vector3d& imgPixel = pixels[i * width + j];

            pixel[0] = imgPixel[2];
            pixel[1] = imgPixel[1];
            pixel[2] = imgPixel[0];
        }
    }

    return mat;
}

cv::Mat Image::toCvMat8UC3() const {
    auto mat = toCvMat();
    cv::Mat mat8uc3;
    mat.convertTo(mat8uc3, CV_8UC3, 255);
    return mat8uc3;
}