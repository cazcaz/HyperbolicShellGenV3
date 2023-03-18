#include "shellGen.h"
#include "circleGen.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include "LBFGS.h"
#include "energyFunction.h"
#include "LineGraphOutputter.h"

using Eigen::Vector3d;
namespace fs = std::filesystem;

ShellGen::ShellGen(ShellParams &parameters) : m_parameters(parameters), m_initLength(parameters.extensionLength), m_radialDist(parameters.radius)
{
    ShellName namer;
    std::string fileName = namer.makeName(m_parameters);
    m_outputDirectory = "./Surfaces/" + fileName;
    fs::create_directory(m_outputDirectory);
};
ShellGen::~ShellGen(){};

void ShellGen::setInitCurve()
{
    Vector3d centre(m_parameters.centreX, m_parameters.centreY, m_parameters.centreZ);
    m_surface = RadialSurface();
    std::vector<Vector3d> initCurve;
    std::vector<Vector3d> secondCurve;
    CircleGen circlemaker;
    circlemaker.makeCircle(m_parameters.radius, centre, m_parameters.resolution, initCurve);
    m_surface.addCurve(initCurve);
}

bool ShellGen::expandCurve()
{
    int curveCount = m_surface.getCurveCount();
    if (curveCount == 0)
    {
        return false;
    }
    std::vector<Vector3d> binormals;
    std::vector<Vector3d> normals;
    std::vector<Vector3d> tangents;
    double initialDist = m_parameters.radius;
    double lon = lengthFunction(m_radialDist, initialDist);
    int nextRingSize = int(lengthFunction(m_radialDist + m_parameters.extensionLength, initialDist) * m_parameters.resolution / (2 * M_PI * initialDist));

    double angleChange = 2 * M_PI / nextRingSize;

    // Keep to remove the need to have an initial expansion before the normal growth
    // Should change nothing
    for (int i = 0; i < nextRingSize; i++)
    {
        double pointParameter = 2 * M_PI * double(i) / double(nextRingSize);
        Vector3d nextTangent, nextNormal;
        if (curveCount == 1)
        {
            nextTangent = Vector3d(-std::sin(pointParameter), std::cos(pointParameter), 0);
            nextNormal = Vector3d(0, 0, 1);
        }
        else
        {
            nextTangent = m_surface.getEdgeTangent(curveCount - 1, pointParameter);
            nextNormal = m_surface.getEdgeNormal(curveCount - 1, pointParameter);
            if (nextNormal[2] == -1)
            {
                nextNormal[2] = 1;
            }
        }
        tangents.push_back(nextTangent);
        normals.push_back(nextNormal);
    }

    for (int i = 0; i < nextRingSize; i++)
    {
        Vector3d nextBinormal(tangents[i][1] * normals[i][2] - tangents[i][2] * normals[i][1], tangents[i][2] * normals[i][0] - tangents[i][0] * normals[i][2], tangents[i][0] * normals[i][1] - tangents[i][1] * normals[i][0]);
        nextBinormal.normalize();
        binormals.push_back(nextBinormal);
    }

    bool success = true;

    // minimsation time
    std::vector<Vector3d> extendedPrevCurve;
    for (int i = 0; i < nextRingSize; i++)
    {
        double pointParameter = 2 * M_PI * double(i) / double(nextRingSize);
        Vector3d nextPoint = m_surface.getPoint(curveCount - 1, pointParameter);
        extendedPrevCurve.push_back(nextPoint);
    }

    EnergyFunction energyFunctional(m_surface, extendedPrevCurve, binormals, normals, m_parameters, m_radialDist + m_parameters.extensionLength, m_outputDirectory);
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = 200;
    LBFGSpp::LBFGSSolver<double> solver(param);
    double energy;
    std::srand((unsigned int) time(0));
    VectorXd input = 0.01 * VectorXd::Random(nextRingSize);
    // for (int ignore = 0; ignore < nextRingSize; ignore++) {
    //     //input[ignore] += 0.01 * std::cos(double(m_parameters.period) * double(ignore)/double(nextRingSize) * M_PI * 2);
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
    try
    {
        int iterCount = solver.minimize(energyFunctional, input, energy);
        m_surface.addIterCount(iterCount);
        if (iterCount == 200)
        {
            // m_parameters.extensionLength *= 0.5;
            std::cout << "Max iterations reached, halving extension length and trying again." << std::endl;
            success = true;
        }
    }
    catch (...)
    {
        std::cout << "Energy increased in minimisation, continuing." << std::endl;
        // return false;
    }

    std::vector<Vector3d> nextCurve;
    // std::vector<Vector3d> testCurve;
    for (int i = 0; i < nextRingSize; i++)
    {
        if (std::isnan(input[i]))
        {
            std::cout << "Failed from nan input." << std::endl;
            return false;
        }
        nextCurve.push_back(extendedPrevCurve[i] + m_parameters.extensionLength * binormals[i] + m_parameters.extensionLength * input[i] * normals[i]);
        // testCurve.push_back(extendedPrevCurve[i]);
    }
    if (success)
    {
        // m_surface.addCurve(testCurve);
        m_surface.addCurve(nextCurve);
        m_recordedExtensionLengths.push_back(m_parameters.extensionLength);
        m_radialDist += m_parameters.extensionLength;
    }
    m_recordedInputs.push_back(input);
    return true;
}

void ShellGen::expandCurveNTimes()
{
    if (m_parameters.expansions == 0)
    {
        int k = 0;
        // set a max to stop it getting ridiculous
        while (expandCurve() && k < 4000)
        {
            k++;
        }
        return;
    }
    else
    {
        for (int iteration = 0; iteration < m_parameters.expansions; iteration++)
        {
            // std::cout << iteration << std::endl;
            // if (iteration % 10) {
            //     printSurface();
            // }
            if (!expandCurve())
            {
                m_parameters.expansions = iteration;
                return;
            }
        }
    }
    m_parameters.expansions = m_surface.getCurveCount();
}

int ShellGen::correctIndex(int index)
{
    if (index >= m_parameters.resolution)
    {
        return correctIndex(index - m_parameters.resolution);
    }
    else if (index < 0)
    {
        return correctIndex(index + m_parameters.resolution);
    }
    return index;
};

void ShellGen::printSurface()
{

    ShellName namer;
    std::string fileName = namer.makeName(m_parameters);
    int surfaceLength = m_surface.getCurveCount();
    if (surfaceLength < 2)
    {
        // not enough information to make a surface
        return;
    }
    std::ofstream surfaceFile(m_outputDirectory + "/surface.txt");
    surfaceFile << m_surface;
    surfaceFile.close();

    // Calculate curvatures of every point to output to m_curvatureDirectory
    if (m_surface.getCurveCount() > 2 && !m_parameters.saveEveryFrame)
    {
        std::ofstream curvatureFile(m_outputDirectory + "/curvature.txt");
        curvatureFile << std::fixed << "Res: " << m_parameters.resolution << std::endl
              << "Exp: " << m_parameters.expansions << std::endl
              << "Len: " << m_parameters.extensionLength << std::endl
              << "SC: " << m_parameters.springCoeff << std::endl
              << "BS: " << m_parameters.bendingStiffness << std::endl
              << "DC:" << m_parameters.desiredCurvature << std::endl << "?";
        int totalVertices = m_surface.surfaceSize();
        double currentMeanCurv;
        double currentGaussCurv;
        for (int currentCurve = 2; currentCurve < m_surface.getCurveCount(); currentCurve++)
        {
            for (int i = m_surface.curveStartIndex(currentCurve - 1); i < m_surface.curveStartIndex(currentCurve); i++)
            {
                m_surface.curvatures(i, currentGaussCurv, currentMeanCurv);
                curvatureFile << currentGaussCurv << " " << currentMeanCurv;
                if (i != m_surface.curveStartIndex(currentCurve) - 1)
                {
                    curvatureFile << ":";
                }
            }
            if (currentCurve != m_surface.getCurveCount() - 1)
            {
                curvatureFile << "|";
            }
        }
        curvatureFile.close();
    }

    // Calculate the lengths of each ring and output it to a length file with the desired length
    if (m_surface.getCurveCount() > 2 && !m_parameters.saveEveryFrame)
    {
        std::ofstream lengthFile(m_outputDirectory + "/lengthProfile.txt");
        std::string fileName = "lengthProfile";
        std::string plotTitle = "Curve Lengths";
        std::string xAxisLabel = "Radius";
        std::string yAxisLabel = "Length";
        std::string currentValuesLegend = "Actual Radius Lengths";
        std::string expValuesLegend = "Desired Radius Lengths";
        int totalCurves = m_surface.getCurveCount();
        double currentLength;
        double expectedLength;
        std::vector<double> currentLengths;
        std::vector<double> expectedLengths;
        double currentRadius = m_parameters.radius;
        std::vector<double> radii;
        for (int i = 0; i < totalCurves; i++)
        {
            currentLengths.push_back(m_surface.getCurveLength(i));
            expectedLengths.push_back(lengthFunction(currentRadius, m_parameters.radius));
            radii.push_back(currentRadius);
            currentRadius += m_recordedExtensionLengths[i];
        }
        LineGraphOutputter lengthCurvePlotter(m_outputDirectory, fileName, plotTitle, xAxisLabel, yAxisLabel, m_parameters);
        lengthCurvePlotter.addXValues(radii);
        lengthCurvePlotter.addData(currentLengths, currentValuesLegend);
        lengthCurvePlotter.addData(expectedLengths, expValuesLegend);
        lengthCurvePlotter.writeData();
    }

    // Output inputs for input graph
    if (m_surface.getCurveCount() > 2 && !m_parameters.saveEveryFrame)
    {
        std::ofstream displacementsFile(m_outputDirectory + "/displacements.txt");
        displacementsFile << std::fixed << "Res: " << m_parameters.resolution << std::endl
              << "Exp: " << m_parameters.expansions << std::endl
              << "Len: " << m_parameters.extensionLength << std::endl
              << "SC: " << m_parameters.springCoeff << std::endl
              << "BS: " << m_parameters.bendingStiffness << std::endl
              << "DC:" << m_parameters.desiredCurvature << std::endl << "?";
        for (int inputsIndex = 0; inputsIndex < m_recordedInputs.size(); inputsIndex++)
        {
            for (int inputIndex = 0; inputIndex < m_recordedInputs[inputsIndex].size(); inputIndex++)
            {
                displacementsFile << std::fixed << m_recordedInputs[inputsIndex][inputIndex];
                if (inputIndex != m_recordedInputs[inputsIndex].size() - 1)
                {
                    displacementsFile << ",";
                }
            }
            if (inputsIndex != m_recordedInputs.size() - 1)
            {
                displacementsFile << "|";
            }
        }
        displacementsFile.close();
    }
}

double ShellGen::lengthFunction(double t, double t0)
{
    double sqrtDC = std::sqrt(std::abs(m_parameters.desiredCurvature));
    return 2 * M_PI * (1 / sqrtDC * std::sinh(sqrtDC * (t - t0)) + t0);
};
