#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"
#include <cstdlib>

using Eigen::Vector3d;
int main(int, char**) {
    system("sudo ./cleanup.sh");
    
    // // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    // std::vector<ShellParams> parameterList;
    // BatchGen massCalcer;
    // double bendStiffChange = (10 - 0.1)/31;
    // double lengthStiffChange = (10 - 0.1)/31;
    // ShellParams parameters;
    // parameters.centreX = 0;
    // parameters.centreY = 0;
    // parameters.centreZ = 0;
    // parameters.desiredCurvature = -0.01;
    // parameters.expansions = 23;
    // parameters.extensionLength = 0.1;


    // for (int k=0; k < 31; k++){
    //     for (int j=0; j < 31; j++)
    //         parameterList.push_back(parameters);
    //         parameters.bendingStiffness += bendStiffChange;
    //         parameters.surfaceIndex += 1;
    //     parameters.lengthStiffness += lengthStiffChange;
    //     parameters.bendingStiffness = 0.1;
    // }

    // // for (int j=0; j < 10; j++){
    // // parameterList.push_back(parameters);
    // //     parameters.period+=1;
    // // }

    // //Push parameters to the parameterList

    // massCalcer.calculateAll(parameterList);

    ShellParams parameters;
    parameters.expansions = 5;
    parameters.desiredCurvature = -1;
    ShellGen shellGenerator(parameters);
    shellGenerator.setInitCurve();
    shellGenerator.expandCurveNTimes();
    shellGenerator.printSurface();

    std::cout << "Surface generation complete, converting to .stl" << std::endl;
    system("sudo python3 txtToStl.py");
    std::cout << "Generating curvature .png s" << std::endl;
    system("sudo python3 curvTxtToPng.py");
    std::cout << "Generating length profile .png s" << std::endl;
    system("sudo python3 lengthTxtToPng.py");
    std::cout << "Generating displacement profile .png s" << std::endl;
    system("sudo python3 dispToPng.py");
    return 0;
}
