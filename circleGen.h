#pragma once

#include <vector>
#include <eigen3/Eigen/Core>

using Eigen::Vector3d;

class CircleGen
{
public:
    CircleGen();
    ~CircleGen();

    void makeCircle(double radius, Vector3d &centre, int resolution, std::vector<Vector3d> &curve);
};