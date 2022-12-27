#include "energyFunction.h"
#include <cmath>
#include <iostream>

using Eigen::Vector3d;
using Eigen::VectorXd;

EnergyFunction::EnergyFunction(RadialSurface& surface,
                               std::vector<Vector3d>& normals,
                               std::vector<Vector3d>& binormals,
                               ShellParams& parameters,
                               double radialDist) :
                               m_surface(surface) , 
                               m_normals(normals) ,
                               m_binormals(binormals) ,
                               m_parameters(parameters) ,
                               m_radialDist(radialDist) {
};

EnergyFunction::~EnergyFunction() = default;

double EnergyFunction::operator()(const VectorXd& inputs, VectorXd& derivatives){
    std::vector<Vector3d> nextCurve;
    int curveCount = m_surface.getCurveCount();
    double sParam;
    int curveSize = m_surface.getCurve(curveCount-1).size();
    for (int i=0; i<curveSize; i++) {
        sParam = 2 * M_PI/curveSize * i;
        nextCurve.push_back(m_surface.getPoint(curveCount-1, sParam) + m_parameters.extensionLength * m_normals[i] + m_parameters.extensionLength * inputs[i] * m_binormals[i]);
    }
    // Here find the mean and gaussian curvatures of the surface at every vertex


    return 0;
};

Vector3d EnergyFunction::normalVecDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db){
    double norm = (a-b).norm();
    return (da-db)/norm - ((da-db).dot(a-b)/std::pow(norm,3)) * (a-b);
};

double EnergyFunction::normDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db){
    double norm = (a-b).norm();
    if (norm == 0) {
        return 0;
    }
    return (a-b).dot(da-db)/norm;
};

double EnergyFunction::dxDij(double x, double xdij, Vector3d p1, Vector3d p2, Vector3d p3, Vector3d dp1, Vector3d dp2, Vector3d dp3){
    double norm1 = (p3-p2).norm();
    double norm2 = (p2-p1).norm();
    return 2/(norm1 + norm2) * (xdij/std::pow(x-1,2))  - (normDeriv(p3,p2,dp3,dp2) + normDeriv(p2,p1,dp2,dp1))/(std::pow(norm1+norm2, 2)) * std::pow(std::tan((M_PI - std::acos(x))/2),2);
};

int EnergyFunction::correctIndex(int index){
    if (index >= m_parameters.resolution) {
        return correctIndex(index - m_parameters.resolution);
    } else if (index < 0) {
        return correctIndex(index + m_parameters.resolution);
    }
    return index;
};

double EnergyFunction::lengthFunction(double t, double t0){
    double sqrtDC = std::sqrt(m_parameters.desiredCurvature);
    return 2 * M_PI * ( 1/sqrtDC * std::sinh(sqrtDC * (t-t0)) + t0);
};
    
double EnergyFunction::rescaleEnergyFunction(double t, double t0){
    return 1/(std::pow(m_parameters.resolution,2) *lengthFunction(t,t0));
};

double EnergyFunction::heavisideApprox(double t){
    return 1/(1+ std::exp(-4*t));
};

double EnergyFunction::inverseLengthFunction(double t, double t0){
    double sqrtDC = std::sqrt(m_parameters.desiredCurvature);
    return t0 + std::asinh(sqrtDC * t/2 * M_PI)/sqrtDC;
};

double EnergyFunction::bendingEnergy(Vector3d a, Vector3d b, Vector3d c){
    double cosAngle = (c-b).normalized().dot((a-b).normalized());
    return (1/((c-b).norm() + (a-b).norm()) * std::pow(std::tan((M_PI - std::acos(cosAngle))/2),2));
};

double EnergyFunction::bendingEnergyDeriv(Vector3d a, Vector3d b, Vector3d c, Vector3d da, Vector3d db, Vector3d dc)
{
    double dxdij = ((c-b).normalized().dot(normalVecDeriv(a,b,da,db)) + (a-b).normalized().dot(normalVecDeriv(c,b,dc,db)));
    double cosAngle = (c-b).normalized().dot((a-b).normalized());
    return dxDij(cosAngle, dxdij, a, b, c, da, db, dc);
};
