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
    // // double meanStiffChange = (0.01 - 0.0001)/10;
    // // double gaussStiffChange = (0.01 - 0.0001)/10;
    // double bendStiffChange = (0.1 - 0.001)/100;
    // ShellParams parameters;
    // parameters.centreX = 0;
    // parameters.centreY = 0;
    // parameters.centreZ = 0;
    // parameters.desiredCurvature = -1;
    // parameters.expansions = 30;
    // parameters.resolution = 500;
    // parameters.extensionLength = 0.1;
    // parameters.bendingStiffness = 0.1;


    //     // for (int k=0; k < 20; k++){
    //     //     parameterList.push_back(parameters);
    //     //     parameters.bendingStiffness += bendStiffChange;
    //     // }

    // for (int j=0; j < 10; j++){
    //     parameterList.push_back(parameters);
    //     parameters.period+=1;
    // }

    // //Push parameters to the parameterList

    // massCalcer.calculateAll(parameterList);

    ShellParams parameters;
    parameters.expansions = 50;
    parameters.desiredCurvature = -1;
    ShellGen shellGenerator(parameters);
    shellGenerator.setInitCurve();
    shellGenerator.expandCurveNTimes();
    shellGenerator.printSurface();


    std::cout << "Surface generation complete, converting to .stl" << std::endl;
    system("sudo python3 txtToStl.py");
    std::cout << "Generating curvature .png s" << std::endl;
    system("sudo python3 curvTxtToPng.py");
    return 0;
}
