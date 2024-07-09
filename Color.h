//
// Created by sty-g on 2023/11/16.
//

#ifndef HELLO_EMBREE_COLOR_H
#define HELLO_EMBREE_COLOR_H

#include <Eigen/Dense>

using Color = Eigen::Vector3d;

double getLuminance(const Color &c);

#endif //HELLO_EMBREE_COLOR_H
