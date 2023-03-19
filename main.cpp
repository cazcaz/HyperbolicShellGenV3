#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"
#include <cstdlib>

using Eigen::Vector3d;
int main(int, char **)
{
    system("sudo ./cleanup.sh");

    // To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    // std::vector<ShellParams> parameterList;
    // BatchGen massCalcer;
    // double DCChange = 0.2;

    // for (int j = 0; j < 10; j++)
    // {
    //     ShellParams parameters;
    //     parameters.expansions = 60;
    //     parameters.desiredCurvature = -1 - double(j)*DCChange;
    //     parameters.surfaceIndex = j;
    //     parameterList.push_back(parameters);
    // }

    // // for (int j=0; j < 10; j++){
    // // parameterList.push_back(parameters);
    // //     parameters.period+=1;
    // // }

    // // Push parameters to the parameterList

    // massCalcer.calculateAll(parameterList);

    ShellParams parameters;
    parameters.expansions = 30;
    parameters.desiredCurvature = -1.5;
    ShellGen shellGenerator(parameters);
    shellGenerator.setInitCurve();
    shellGenerator.expandCurveNTimes();
    shellGenerator.printSurface();

    std::cout << "Surface generation complete, converting to .stl" << std::endl;
    system("sudo python3 txtToStl.py");
    std::cout << "Generating polar graph .png s" << std::endl;
    system("sudo python3 polarGraph.py");
    std::cout << "Generating line graphs .png s" << std::endl;
    system("sudo python3 lineGraphPng.py");
    return 0;
}
