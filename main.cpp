#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"

using Eigen::Vector3d;
int main(int, char**) {
    
    // // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    // std::vector<ShellParams> parameterList;
    // BatchGen massCalcer;
    // double stiffChange = (10000-1)/10;
    // double curveChange = (0.0001 - 0.000001)/10;
    // double stiffness = 1;

    // ShellParams parameters;
    // parameters.centreX = 0;
    // parameters.centreY = 0;
    // parameters.centreZ = 0;
    // parameters.desiredCurvature = 0.000001;
    // parameters.expansions = 15;
    // parameters.resolution = 100;
    // parameters.extensionLength = 0.1;
    // parameters.stiffnessRatio = 100;


    // for (int j=0; j < 10; j++){
    //     for (int k=0; k < 10; k++){
    //         parameterList.push_back(parameters);
    //         parameters.desiredCurvature += curveChange;
    //     }
    //     parameters.desiredCurvature = 0.000001;
    // parameters.stiffnessRatio += stiffChange;
    // }

    // //Push parameters to the parameterList

    // massCalcer.calculateAll(parameterList);

    ShellParams parameters;
    parameters.expansions = 50;
    parameters.desiredCurvature = -0.001;
    ShellGen shellGenerator(parameters);
    shellGenerator.setInitCurve();
    shellGenerator.expandCurveNTimes();
    shellGenerator.printSurface();

    return 0;
}
