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
        for (int j=0; j < 101; j++){
            desiredCurv = 0.5-0.01*j;
            for (int i=0; i < 11; i++) {
                ShellParams parameters;
                parameters.desiredCurvature = desiredCurv;
                parameters.strainCoeff = 10*i;
                parameterList.push_back(parameters);
            }
        }

    massCalcer.calculateAll(parameterList);

    // ShellParams parameters;
    // parameters.expansions = 10;
    // parameters.desiredCurvature = -0.1;
    // ShellGen shellGenerator(parameters);
    // shellGenerator.setInitCurve();
    // shellGenerator.expandCurveNTimes();
    // shellGenerator.printSurface();

    return 0;
}
