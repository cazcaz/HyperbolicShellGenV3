#include "shellGen.h"
#include "circleGen.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "LBFGS.h"
#include "energyFunction.h"

using Eigen::Vector3d;
namespace fs = std::filesystem;

ShellGen::ShellGen(ShellParams& parameters) : m_parameters(parameters) , m_initLength(parameters.extensionLength), m_radialDist(parameters.radius + m_parameters.extensionLength){
    ShellName namer;
    std::string fileName = namer.makeName(m_parameters);
    m_outputDirectory ="./Surfaces/" + fileName;
    fs::create_directory(m_outputDirectory);
};
ShellGen::~ShellGen() {};

void ShellGen::setInitCurve() {
    Vector3d centre(m_parameters.centreX, m_parameters.centreY, m_parameters.centreZ);
    m_surface = RadialSurface();
    std::vector<Vector3d> initCurve;
    std::vector<Vector3d> secondCurve;
    CircleGen circlemaker;
    circlemaker.makeCircle(m_parameters.radius, centre, m_parameters.resolution, initCurve);
    int nextRingRes = lengthFunction(m_parameters.radius + m_initLength, m_parameters.radius) * m_parameters.resolution / (2 * M_PI * m_parameters.radius);
    circlemaker.makeCircle(m_parameters.radius + m_initLength, centre, nextRingRes, secondCurve);
    m_surface.addCurve(initCurve);
    m_surface.addCurve(secondCurve);
    m_recordedExtensionLengths.push_back(m_parameters.extensionLength);
}

bool ShellGen::expandCurve() {
    int curveCount = m_surface.getCurveCount();
    if(curveCount == 0) {
        return false;
    }
    std::vector<Vector3d> binormals;
    std::vector<Vector3d> normals;
    std::vector<Vector3d> tangents;
    double initialDist = m_parameters.radius;
    double lon = lengthFunction(m_radialDist, initialDist);
    int nextRingSize = int(lengthFunction(m_radialDist + m_parameters.extensionLength, initialDist) * m_parameters.resolution / (2 * M_PI * initialDist));

    double angleChange = 2 * M_PI / nextRingSize;
    for (int i =0; i<nextRingSize;i++){
            double pointParameter = 2 * M_PI * double(i)/double(nextRingSize);
            Vector3d nextTangent = m_surface.getEdgeTangent(curveCount-1, pointParameter);
            Vector3d nextNormal = m_surface.getEdgeNormal(curveCount-1, pointParameter);
            if (nextNormal[2] == -1) {
                nextNormal[2] = 1;
            }
            tangents.push_back(nextTangent);
            normals.push_back(nextNormal);
    }
    
    if (curveCount == 1) {
        for (Vector3d firstCurvePoint : m_surface.getCurve(1)){
            binormals.push_back(firstCurvePoint.normalized());
        }
    } else {
        for (int i =0; i<nextRingSize;i++){
            Vector3d nextBinormal(tangents[i][1]*normals[i][2] - tangents[i][2]*normals[i][1], tangents[i][2]*normals[i][0] - tangents[i][0]*normals[i][2], tangents[i][0]*normals[i][1] - tangents[i][1]*normals[i][0]);
            nextBinormal.normalize();
            binormals.push_back(nextBinormal);
        }
    }
    //Nicely behaved tangents
    

    bool success = true;

    //minimsation time
    std::vector<Vector3d> extendedPrevCurve;
    for (int i = 0; i < nextRingSize; i++) {
        double pointParameter = 2 * M_PI * double(i)/double(nextRingSize);
        Vector3d nextPoint = m_surface.getPoint(curveCount-1, pointParameter);
        extendedPrevCurve.push_back(nextPoint);
    }

    EnergyFunction energyFunctional(m_surface, extendedPrevCurve, binormals, normals, m_parameters, m_radialDist + m_parameters.extensionLength, m_outputDirectory);
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = 200;
    LBFGSpp::LBFGSSolver<double> solver(param);
    double energy;
    VectorXd input = 0.5 * VectorXd::Random(nextRingSize);
    // for (int ignore = 0; ignore < nextRingSize; ignore++) {
    //     input[ignore] += 0.05 * std::cos(double(m_parameters.period) * double(ignore)/double(nextRingSize) * M_PI * 2);
    // }

    // Used to make a linear approx. of the derivative for testing
    // VectorXd inputChanged = input;
    // double h = 0.00000001;
    // inputChanged[10] += h;
    // VectorXd derivatives = VectorXd::Zero(nextRingSize);
    // double energy2;
    // energy2 = energyFunctional(inputChanged, derivatives);
    // energy = energyFunctional(input, derivatives);
    // std::cout << "Approx: " << (energy2 - energy)/h<< std::endl;
    // std::cout << "Real: " << derivatives[10] << std::endl;
    // double tempEnergyComp = ((energy2 - energy)/h)/derivatives[10];
    // std::cout << "Resolution: " << m_parameters.resolution << " NextRingSize: " << nextRingSize << " Energy Ratio: " <<((energy2 - energy)/h)/derivatives[10] << std::endl;
    
    // std::cout << derivatives.transpose() << std::endl;
    try {
        int iterCount = solver.minimize(energyFunctional, input, energy);
        m_surface.addIterCount(iterCount);
        if (iterCount == 200) {
            //m_parameters.extensionLength *= 0.5;
            std::cout << "Max iterations reached, halving extension length and trying again." << std::endl;
            success = true;
        }
    } catch(...) {
        std::cout << "Failed from error in calculation." << std::endl;
        return false;
    }

    std::vector<Vector3d> nextCurve;
    //std::vector<Vector3d> testCurve;
    for (int i =0; i<nextRingSize;i++) {
        if (std::isnan(input[i])){
            std::cout << "Failed from nan input." << std::endl;
            return false;
        }
        nextCurve.push_back(extendedPrevCurve[i] + m_parameters.extensionLength * binormals[i] + m_parameters.extensionLength * input[i] * normals[i]);
        //testCurve.push_back(extendedPrevCurve[i]);
    }
    if (success) {
        //m_surface.addCurve(testCurve);
        m_surface.addCurve(nextCurve);
        m_recordedExtensionLengths.push_back(m_parameters.extensionLength);
        m_radialDist += m_parameters.extensionLength;
    }
    m_recordedInputs.push_back(input);
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
            // std::cout << iteration << std::endl;
            // if (iteration % 10) {
            //     printSurface();
            // }
            if (!expandCurve()){
                m_parameters.expansions = iteration;
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
    std::ofstream surfaceFile(m_outputDirectory +"/surface.txt");
    surfaceFile << m_surface;
    surfaceFile.close();


    //Calculate curvatures of every point to output to m_curvatureDirectory
    if (m_surface.getCurveCount() > 2 && !m_parameters.saveEveryFrame){
    std::ofstream curvatureFile(m_outputDirectory +"/curvature.txt");
    int totalVertices = m_surface.surfaceSize();
    double currentMeanCurv;
    double currentGaussCurv;
    for (int currentCurve = 2; currentCurve < m_surface.getCurveCount(); currentCurve++){
        for (int i = m_surface.curveStartIndex(currentCurve-1); i < m_surface.curveStartIndex(currentCurve); i++){
            m_surface.curvatures(i, currentGaussCurv, currentMeanCurv);
            curvatureFile << currentGaussCurv << " " << currentMeanCurv;
            if (i != m_surface.curveStartIndex(currentCurve) - 1){
            curvatureFile << ":";
            }
        }
        if (currentCurve != m_surface.getCurveCount() - 1){
            curvatureFile << "|";
        }
    }
    curvatureFile.close();
    }

    //Calculate the lengths of each ring and output it to a length file with the desired length
    if (m_surface.getCurveCount() > 2 && !m_parameters.saveEveryFrame){
    std::ofstream lengthFile(m_outputDirectory +"/lengthProfile.txt");
    int totalCurves = m_surface.getCurveCount();
    double currentLength;
    double expectedLength;
    double currentRadius = m_parameters.radius;
    for (int i =  0; i < totalCurves; i++){
        currentLength = m_surface.getCurveLength(i);
        expectedLength = lengthFunction(currentRadius,m_parameters.radius);
        currentRadius += m_recordedExtensionLengths[i];
        lengthFile << currentRadius << " " << currentLength << " " << expectedLength;
        if (i != totalCurves-1){
        lengthFile << ":";
        }
    }
    lengthFile.close();
    }

    //Output inputs for input graph
    if (m_surface.getCurveCount() > 2 && !m_parameters.saveEveryFrame){
    std::ofstream displacementsFile(m_outputDirectory +"/displacements.txt");
    for (int inputsIndex = 0; inputsIndex<m_recordedInputs.size(); inputsIndex++) {
        for (int inputIndex=0; inputIndex < m_recordedInputs[inputsIndex].size();inputIndex++) {
            displacementsFile << std::fixed << m_recordedInputs[inputsIndex][inputIndex];
            if (inputIndex != m_recordedInputs[inputsIndex].size()-1) {
                displacementsFile << ",";
            }
        }
        if (inputsIndex != m_recordedInputs.size()-1) {
            displacementsFile << "|";
        }
    }
    displacementsFile.close();
    }
}

double ShellGen::lengthFunction(double t, double t0){
    double sqrtDC = std::sqrt(std::abs(m_parameters.desiredCurvature));
    return 2 * M_PI * ( 1/sqrtDC * std::sinh(sqrtDC * (t-t0)) + t0);
};
