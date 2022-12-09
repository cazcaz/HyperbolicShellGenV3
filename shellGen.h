#pragma once
#include <memory>
#include <vector>
#include <Eigen/core>
#include <string>
#include "shellParams.h"
#include "shellNaming.h"
using Eigen::Vector3d;
using Eigen::VectorXd;

class ShellGen {
    public:
        ShellGen(ShellParams& parameters);
        ~ShellGen();

        void setInitCurve();

        void expandCurve();

        void expandCurveNTimes(int iterations);

        void printSurface();

    private:
        ShellParams& m_parameters;
        std::vector<std::vector<Vector3d>> m_surface;
        VectorXd m_prevSol;
};