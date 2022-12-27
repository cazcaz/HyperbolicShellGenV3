#include "batchGen.h"
#include "surface.h"
#include "circleGen.h"

using Eigen::Vector3d;
int main(int, char**) {
    
    //To use, populate parameterList with parameters of surfaces to be calculated, then call batchGen to calcuate them

    // std::vector<ShellParams> parameterList;
    // BatchGen massCalcer;
    
    // //Push parameters to the parameterList
    // ShellParams parameters1;
    // parameters1.expansions = 100;
    // parameterList.push_back(parameters1);

    // massCalcer.calculateAll(parameterList);

    CircleGen circleMaker;
    Surface testSurface;
    Vector3d centre = Vector3d(0,0,1);
    testSurface.addPoint(centre);
    std::vector<Vector3d> circle;
    circleMaker.makeCircle(1, Vector3d::Zero(), 5, circle); 
    for (Vector3d point : circle) {
        testSurface.addPoint(point);
    }
    testSurface.addTriangle(Triangle(0,1,2));
    testSurface.addTriangle(Triangle(0,2,3));
    testSurface.addTriangle(Triangle(0,3,4));
    testSurface.addTriangle(Triangle(0,4,5));
    testSurface.addTriangle(Triangle(0,5,1));
    std::cout << testSurface.meanCurvature(0) << std::endl;
    return 0;
}
