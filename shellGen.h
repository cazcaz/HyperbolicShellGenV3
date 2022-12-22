#pragma once
#include <memory>
#include <vector>
#include <Eigen/core>
#include <string>
#include "shellParams.h"
#include "shellNaming.h"
#include "radialSurface.h"
using Eigen::Vector3d;
using Eigen::VectorXd;

class ShellGen {
    public:
        ShellGen(ShellParams& parameters);
        ~ShellGen();

        void setInitCurve();

        bool expandCurve();

        void expandCurveNTimes();

        int correctIndex(int index);

        void printSurface();

    private:
        struct ShellParams& m_parameters;
        RadialSurface m_surface;
};