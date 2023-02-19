#pragma once
#include <memory>
#include <vector>
#include <eigen3/Eigen/Core>
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

        double lengthFunction(double t, double t0);

        void printSurface();

    private:
        struct ShellParams& m_parameters;
        RadialSurface m_surface;
        double m_initLength;
        double m_radialDist;
        std::string m_outputDirectory = "";
};