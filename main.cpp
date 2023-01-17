#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"

using Eigen::Vector3d;
int main(int, char**) {
    
    // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    std::vector<ShellParams> parameterList;
    BatchGen massCalcer;
    double desiredCurv;
    //Push parameters to the parameterList
    for (int i=0; i < 26; i++){
        desiredCurv = -0.0001 - 0.1 * i;
        for (int j=0; j < 26; j++){
            ShellParams parameters;
            parameters.desiredCurvature = desiredCurv;
            parameters.stiffnessRatio = 0.001 + 0.1 * j;
            parameterList.push_back(parameters);
        }

    }

    massCalcer.calculateAll(parameterList);

    // ShellParams parameters;
    // parameters.expansions = 10;
    // parameters.desiredCurvature = 0;
    // ShellGen shellGenerator(parameters);
    // shellGenerator.setInitCurve();
    // shellGenerator.expandCurveNTimes();
    // shellGenerator.printSurface();

    return 0;
}
