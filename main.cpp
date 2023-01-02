#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"

using Eigen::Vector3d;
int main(int, char**) {
    
    // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    std::vector<ShellParams> parameterList;
    BatchGen massCalcer;
    
    //Push parameters to the parameterList
    for (int i=0; i < 100; i++){
        ShellParams parameters;
        parameters.expansions = 20;
        parameters.desiredCurvature = 0.0005 - 0.00001 * i;
        parameterList.push_back(parameters);
    }

    massCalcer.calculateAll(parameterList);

    return 0;
}
