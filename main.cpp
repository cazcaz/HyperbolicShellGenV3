#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"
#include <cstdlib>

using Eigen::Vector3d;
int main(int, char**) {
    
    // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    std::vector<ShellParams> parameterList;
    BatchGen massCalcer;
    // double meanStiffChange = (0.01 - 0.0001)/10;
    // double gaussStiffChange = (0.01 - 0.0001)/10;
    double bendStiffChange = (0.1 - 0.001)/100;

    ShellParams parameters;
    parameters.centreX = 0;
    parameters.centreY = 0;
    parameters.centreZ = 0;
    parameters.desiredCurvature = -1;
    parameters.expansions = 4;
    parameters.resolution = 100;
    parameters.extensionLength = 0.1;
    // parameters.meanStiffness = 0.0001;
    // parameters.gaussStiffness = 0.0001;
    parameters.bendingStiffness = 0.001;


    // for (int j=0; j < 10; j++){
    //     for (int k=0; k < 10; k++){
    //         parameterList.push_back(parameters);
    //         parameters.gaussStiffness += gaussStiffChange;
    //     }
    //     parameters.gaussStiffness = 0.0001;
    // parameters.meanStiffness += meanStiffChange;
    // }

    for (int j=0; j < 3; j++){
        parameterList.push_back(parameters);
        parameters.bendingStiffness += bendStiffChange;
    }

    //Push parameters to the parameterList

    massCalcer.calculateAll(parameterList);

    // ShellParams parameters;
    // parameters.expansions = 3;
    // parameters.desiredCurvature = -1;
    // ShellGen shellGenerator(parameters);
    // shellGenerator.setInitCurve();
    // shellGenerator.expandCurveNTimes();
    // shellGenerator.printSurface();

    system("sudo python3 txtToStl.py");
    return 0;
}
