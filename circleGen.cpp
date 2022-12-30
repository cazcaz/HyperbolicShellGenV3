#pragma once
#include "circleGen.h"
#include <cmath>
#include <iostream>


using  Eigen::Vector3d;

CircleGen::CircleGen() {
    
};

CircleGen::~CircleGen() {};

void CircleGen::makeCircle(double radius, Vector3d& centre, int resolution, std::vector<Vector3d> &curve) {
    for (int i = 0;  i < resolution; i++){
        float angle = 2 * M_PI * i / resolution;
        Vector3d next(radius * std::cos(angle), radius * std::sin(angle), 0);
        curve.push_back(centre + next);
    }
}