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
        Vector3d normalVecDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        double normDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        Vector3d normVecDeriv(Vector3d& a, Vector3d& da);
        double atan2Deriv(double x, double y, double dx, double dy);
        Vector3d crossProd(Vector3d &a, Vector3d &b);
        double angleDeriv(double angle, Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &da, Vector3d &db, Vector3d &dc);
        double dotDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        Vector3d crossDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        void dihedralAngleDeriv(Vector3d& a, Vector3d& b, Vector3d& c, Vector3d& d, Vector3d& da, Vector3d& db, Vector3d& dc, Vector3d& dd, double& angleResult, double& derivResult);
        double triangleAreaDeriv(double angle1, double angle2, double dangle1, double dangle2 ,Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &d, Vector3d &da, Vector3d &db, Vector3d &dc, Vector3d& dd, double& dl1, double& dl2, double& dl3);
        double lengthFunction(double t, double t0);
        double bendingEnergy(Vector3d a, Vector3d b, Vector3d c);
        double bendingEnergyDeriv(Vector3d a, Vector3d b, Vector3d c, Vector3d da, Vector3d db, Vector3d dc);
        double dxDij(double x, double xdij, Vector3d p1, Vector3d p2, Vector3d p3, Vector3d dp1, Vector3d dp2, Vector3d dp3);
        double rescaleEnergyFunction(double t, double t0, int nextCurveLength);
    private:
        RadialSurface m_surface;
        std::vector<Vector3d> m_normals;
        std::vector<Vector3d> m_binormals;
        std::vector<Vector3d> m_prevCurve;
        struct ShellParams& m_parameters;
        double m_radialDist;
        bool m_firstRun;
        std::vector<double> m_origTriangleSizes;
        std::vector<std::vector<double>> m_origTrianglePenaltySizes;
};