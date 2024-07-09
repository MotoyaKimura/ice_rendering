//
// Created by sty-g on 2023/11/16.
//

#include "Color.h"

double getLuminance(const Color &c) {
    return c.dot(Eigen::Vector3d(0.2126, 0.7152, 0.0722));
}