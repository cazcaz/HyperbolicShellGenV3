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
    parameters.expansions = 3;
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
