#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"

using Eigen::Vector3d;
int main(int, char**) {
    
    // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    std::vector<ShellParams> parameterList;
    BatchGen massCalcer;
    
    //Push parameters to the parameterList
    for (int i=0; i < 1; i++){
        ShellParams parameters;
        parameters.expansions = 20;
        parameters.desiredCurvature = -1 - 0.01 * i;
        parameterList.push_back(parameters);
    }

    massCalcer.calculateAll(parameterList);

    // ShellParams parameters;
    // parameters.expansions = 20;
    // parameters.desiredCurvature = 0.0005;
    // ShellGen shellGenerator(parameters);
    // shellGenerator.setInitCurve();
    // shellGenerator.expandCurveNTimes();

    return 0;
}
