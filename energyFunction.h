#pragma once

#include <eigen3/Eigen/Core>
#include <vector>
#include "shellParams.h"
#include "radialSurface.h"

using Eigen::Vector3d;
using Eigen::VectorXd;

class EnergyFunction {
    public:
        EnergyFunction(RadialSurface& surface, std::vector<Vector3d>& extendedPrevCurve, std::vector<Vector3d>& normals, std::vector<Vector3d>& binormals, ShellParams& parameters, double radialDist);
        ~EnergyFunction();

        double operator()(const VectorXd& inputs, VectorXd& derivatives);
        double evalEnergy(RadialSurface& extendedSurface);
    private:
        RadialSurface m_surface;
        std::vector<Vector3d> m_normals;
        std::vector<Vector3d> m_binormals;
        std::vector<Vector3d> m_prevCurve;
        struct ShellParams& m_parameters;
        double m_radialDist;
};