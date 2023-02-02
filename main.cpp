#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"

using Eigen::Vector3d;
int main(int, char**) {
    
    // // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    // std::vector<ShellParams> parameterList;
    // BatchGen massCalcer;
    // double meanStiffChange = (0.01 - 0.0001)/10;
    // double gaussStiffChange = (0.01 - 0.0001)/10;

    // ShellParams parameters;
    // parameters.centreX = 0;
    // parameters.centreY = 0;
    // parameters.centreZ = 0;
    // parameters.desiredCurvature = -0.01;
    // parameters.expansions = 1000;
    // parameters.resolution = 30;
    // parameters.extensionLength = 0.1;
    // parameters.meanStiffness = 0.0001;
    // parameters.gaussStiffness = 0.0001;


    // for (int j=0; j < 10; j++){
    //     for (int k=0; k < 10; k++){
    //         parameterList.push_back(parameters);
    //         parameters.gaussStiffness += gaussStiffChange;
    //     }
    //     parameters.gaussStiffness = 0.0001;
    // parameters.meanStiffness += meanStiffChange;
    // }

    // //Push parameters to the parameterList

    // massCalcer.calculateAll(parameterList);

    ShellParams parameters;
    parameters.expansions = 100;
    parameters.desiredCurvature = -0.01;
    ShellGen shellGenerator(parameters);
    shellGenerator.setInitCurve();
    shellGenerator.expandCurveNTimes();
    shellGenerator.printSurface();

    return 0;
}
