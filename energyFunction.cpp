#include "energyFunction.h"
#include <cmath>
#include <iostream>

using Eigen::Vector3d;
using Eigen::VectorXd;

EnergyFunction::EnergyFunction(RadialSurface& surface,
                               std::vector<Vector3d>& extendedPrevCurve,
                               std::vector<Vector3d>& normals,
                               std::vector<Vector3d>& binormals,
                               ShellParams& parameters,
                               double radialDist) :
                               m_surface(surface) , 
                               m_prevCurve(extendedPrevCurve) ,
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
    int curveSize = inputs.size();
    for (int i=0; i<curveSize; i++) {
        nextCurve.push_back(m_prevCurve[i] + m_parameters.extensionLength * m_normals[i] + m_parameters.extensionLength * inputs[i] * m_binormals[i]);
    }
    double totalEnergy = evalEnergy(m_surface);
    double h = 1e-300;
    RadialSurface nextSurface(m_surface);
    nextSurface.addCurve(nextCurve);
    for (int i=0; i<curveSize; i++) {
        Vector3d originalCurvePoint = nextCurve[i];
        Vector3d changedCurvePoint = originalCurvePoint + h * m_parameters.extensionLength * m_binormals[i];
        double beforeEnergy = evalEnergy(nextSurface);
        nextSurface.changeCurvePoint(curveCount, i, changedCurvePoint);
        derivatives[i] = (evalEnergy(nextSurface) - beforeEnergy)/h;
        nextSurface.changeCurvePoint(curveCount, i, originalCurvePoint);
    }
    return totalEnergy;
};

double EnergyFunction::evalEnergy(RadialSurface& extendedSurface)
{
    int curveCount = extendedSurface.getCurveCount();
    // Here find the mean and gaussian curvatures of the surface at every vertex
    double totalGauss = 0;
    double totalMean = 0;
    for (int i = 0; i < extendedSurface.getCurve(curveCount-2).size(); i++) {
        int index = extendedSurface.curveStartIndex(curveCount-2) + i;
        totalGauss += extendedSurface.gaussCurvature(index);
        totalMean += m_parameters.stiffnessRatio / 2 * std::pow(extendedSurface.gaussCurvature(index) - m_parameters.desiredCurvature,2);
    }
    // std::cout << "Mean: " << totalMean << std::endl;
    // std::cout << "Gauss: " << totalGauss << std::endl;
    return totalMean;// + totalGauss;
};

