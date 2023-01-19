#include "shellGen.h"
#include "circleGen.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "LBFGS.h"
#include "energyFunction.h"

using Eigen::Vector3d;

ShellGen::ShellGen(ShellParams& parameters) : m_parameters(parameters) {};
ShellGen::~ShellGen() {};

void ShellGen::setInitCurve() {
    Vector3d centre(m_parameters.centreX, m_parameters.centreY, m_parameters.centreZ);
    m_surface = RadialSurface();
    std::vector<Vector3d> initCurve;
    std::vector<Vector3d> secondCurve;
    CircleGen circlemaker;
    circlemaker.makeCircle(1, centre, m_parameters.resolution, initCurve);
    circlemaker.makeCircle(1 + m_parameters.extensionLength, centre, m_parameters.resolution, secondCurve);
    m_surface.addCurve(initCurve);
    m_surface.addCurve(secondCurve);
}

bool ShellGen::expandCurve() {
    int curveCount = m_surface.getCurveCount();
    if(curveCount == 0) {
        return false;
    }
    std::vector<Vector3d> normals;
    std::vector<Vector3d> binormals;
    std::vector<Vector3d> tangents;
    double initialDist = 1;
    double radialDist = 1 + (curveCount-1) * m_parameters.extensionLength;
    int nextRingSize = int(M_1_PI/(std::asin(1/radialDist * std::sin(M_1_PI/m_parameters.resolution))))/4;
    if (nextRingSize < m_parameters.resolution) {
        nextRingSize = m_parameters.resolution;
    }
    if (curveCount == 1) {
        for (Vector3d firstCurvePoint : m_surface.getCurve(0)){
            normals.push_back(firstCurvePoint.normalized());
        }
    } else {
        for (int i =0; i<nextRingSize;i++){
            double pointParameter = 2 * M_PI * double(i)/double(nextRingSize);
            //Vector3d nextPoint = m_surface.getPoint(curveCount-1, pointParameter) - m_surface.getPoint(curveCount-2, pointParameter);
            Vector3d nextPoint(std::cos(pointParameter), std::sin(pointParameter), 0);
            nextPoint.normalize();
            normals.push_back(nextPoint);
        }
    }
    //Nicely behaved tangents
    double angleChange = 2 * M_PI / nextRingSize;
    for (int i =0; i<nextRingSize;i++){
            Vector3d nextTangent(-std::sin(angleChange * i), std::cos(angleChange *i), 0);
            nextTangent.normalize();
            Vector3d nextBinormal(normals[i][1]*nextTangent[2] - normals[i][2]*nextTangent[1] ,normals[i][2]*nextTangent[0] - normals[i][0]*nextTangent[2], normals[i][0]*nextTangent[1] - normals[i][1]*nextTangent[0]);
            nextBinormal.normalize();
            if (nextBinormal[2] == -1) {
                nextBinormal[2] = 1;
            }
            tangents.push_back(nextTangent);
            binormals.push_back(nextBinormal);
    }

    bool success = true;

    //minimsation time
    std::vector<Vector3d> extendedPrevCurve;
    for (int i = 0; i < nextRingSize; i++) {
        double pointParameter = 2 * M_PI * double(i)/double(nextRingSize);
        Vector3d nextPoint = m_surface.getPoint(curveCount-1, pointParameter);
        extendedPrevCurve.push_back(nextPoint);
    }
    EnergyFunction energyFunctional(m_surface, extendedPrevCurve, normals, binormals, m_parameters, radialDist);
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = 200;
    LBFGSpp::LBFGSSolver<double> solver(param);
    double energy;
    VectorXd input = 0.005 * VectorXd::Random(nextRingSize);

    // // Used to make a linear approx. of the derivative for testing
    // VectorXd inputChanged = input;
    // double h = 0.00001;
    // inputChanged[10] += h;
    // VectorXd derivatives = VectorXd::Zero(nextRingSize);
    // double energy2;
    // energy2 = energyFunctional(inputChanged, derivatives);
    // energy = energyFunctional(input, derivatives);
    // std::cout << "Approx: " << (energy2 - energy)/h<< std::endl;
    // std::cout << "Real: " << derivatives[10] << std::endl;
    // std::cout << derivatives.transpose() << std::endl;
    try {
        int iterCount = solver.minimize(energyFunctional, input, energy);
        if (iterCount == 100) {
            //m_parameters.extensionLength *= 0.5;
            //std::cout << "Max iterations reached, halving extension length and trying again." << std::endl;
            //success = false;
        }
    } catch(...) {
        //std::cout << "Failed from error in calculation." << std::endl;
        return false;
    }

    std::vector<Vector3d> nextCurve;
    for (int i =0; i<nextRingSize;i++) {
        if (std::isnan(input[i])){
            //std::cout << "Failed from nan input." << std::endl;
            return false;
        }
        nextCurve.push_back(extendedPrevCurve[i] + m_parameters.extensionLength * normals[i] + m_parameters.extensionLength * input[i] * binormals[i]);
    }
    if (success) {
        m_surface.addCurve(nextCurve);
    }
    return true;
}

void ShellGen::expandCurveNTimes() {
    if (m_parameters.expansions == 0) {
        int k = 0;
        //set a max to stop it getting ridiculous
        while (expandCurve() && k < 4000) {
            k++;
        }
        return;
    } else {
        for (int iteration = 0; iteration < m_parameters.expansions; iteration++){
            if (!expandCurve()){
                m_parameters.expansions = iteration+1;
                return;
            }
        }
    }
    m_parameters.expansions = m_surface.getCurveCount();
}

int ShellGen::correctIndex(int index){
    if (index >= m_parameters.resolution) {
        return correctIndex(index - m_parameters.resolution);
    } else if (index < 0) {
        return correctIndex(index + m_parameters.resolution);
    }
    return index;
};

void ShellGen::printSurface() {
    
    ShellName namer;
    std::string fileName = namer.makeName(m_parameters);
    int surfaceLength = m_surface.getCurveCount();
    if (surfaceLength < 2) {
        //not enough information to make a surface
        return; 
    }
    std::string path = "/../OutputSurfaceTxts/";
    std::ofstream open(path);
    std::ofstream surfaceFile("../OutputSurfaceTxts/" + fileName + ".txt");
    surfaceFile << m_surface;
    surfaceFile.close();
}


