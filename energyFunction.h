#pragma once

#include <eigen3/Eigen/Core>
#include <vector>
#include "shellParams.h"
#include "radialSurface.h"
#include <ostream>
#include <fstream>

using Eigen::Vector3d;
using Eigen::VectorXd;

class EnergyFunction {
    public:
        EnergyFunction(RadialSurface& surface, std::vector<Vector3d>& extendedPrevCurve, std::vector<Vector3d>& normals, std::vector<Vector3d>& binormals, ShellParams& parameters, double radialDist, std::string outputDirectory);
        ~EnergyFunction();

        double operator()(const VectorXd& inputs, VectorXd& derivatives);
        Vector3d normalVecDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        double normDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        Vector3d normVecDeriv(Vector3d& a, Vector3d& da);
        Vector3d crossProd(Vector3d &a, Vector3d &b);
        double dotDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        Vector3d crossDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db);
        double lengthFunction(double t, double t0);
        double bendingEnergy(Vector3d a, Vector3d b, Vector3d c);
        double bendingEnergyDeriv(Vector3d a, Vector3d b, Vector3d c, Vector3d da, Vector3d db, Vector3d dc);
        double rescaleEnergyFunction(double t, double t0, int nextCurveLength);

        //Currently unused functions kept out of the way to avoid clutter
        double atan2Deriv(double x, double y, double dx, double dy);
        double angleDeriv(double angle, Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &da, Vector3d &db, Vector3d &dc);
        void dihedralAngleDeriv(Vector3d& a, Vector3d& b, Vector3d& c, Vector3d& d, Vector3d& da, Vector3d& db, Vector3d& dc, Vector3d& dd, double& angleResult, double& derivResult);
        double triangleAreaDeriv(double angle1, double angle2, double dangle1, double dangle2 ,Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &d, Vector3d &da, Vector3d &db, Vector3d &dc, Vector3d& dd, double& dl1, double& dl2, double& dl3);

    private:
        RadialSurface m_surface;
        std::vector<Vector3d> m_binormals;
        std::vector<Vector3d> m_normals;
        std::vector<Vector3d> m_prevCurve;
        struct ShellParams& m_parameters;
        double m_radialDist;
        bool m_firstRun;
        std::vector<double> m_origTriangleSizes;
        std::vector<std::vector<double>> m_origTrianglePenaltySizes;
        std::string m_outDirectory;
};