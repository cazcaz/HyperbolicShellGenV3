#include <iostream>
#include "shellGen.h"
#include <Eigen/core>

using Eigen::Vector3d;

int main(int, char**) {
    ShellGen shellGenerator;
    shellGenerator.setInitCurve(1, Vector3d(0,0,0), 5);
    shellGenerator.expandCurve(1, 0.05, 100);
    //shellGenerator.expandCurveNTimes(5,1, 0.05, 100);
    return 0;
}
