#pragma once

#include <eigen3/Eigen/Core>
#include <vector>
#include "shellParams.h"
#include "radialSurface.h"

using Eigen::Vector3d;
using Eigen::VectorXd;

class EnergyFunction {
    public:
        EnergyFunction(RadialSurface& surface, std::vector<Vector3d>& normals, std::vector<Vector3d>& binormals, ShellParams& parameters, double radialDist);
        ~EnergyFunction();

        double operator()(const VectorXd& inputs, VectorXd& derivatives);
        Vector3d normalVecDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        double normDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        double dxDij(double x, double xdij, Vector3d p1, Vector3d p2, Vector3d p3, Vector3d dp1, Vector3d dp2, Vector3d dp3);
        double evalEnergy(RadialSurface& extendedSurface);
        int correctIndex(int index);
        double lengthFunction(double t, double t0);
        double rescaleEnergyFunction(double t, double t0);
        double heavisideApprox(double t);
        double inverseLengthFunction(double t, double t0);
        double bendingEnergy(Vector3d a, Vector3d b, Vector3d c);
        double bendingEnergyDeriv(Vector3d a, Vector3d b, Vector3d c, Vector3d da, Vector3d db, Vector3d dc);
        double deltaFunc(int i, int j);
    private:
        RadialSurface m_surface;
        std::vector<Vector3d> m_normals;
        std::vector<Vector3d> m_binormals;
        struct ShellParams& m_parameters;
        double m_radialDist;
};