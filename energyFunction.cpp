#include "energyFunction.h"
#include <cmath>
#include <iostream>

using Eigen::Vector3d;
using Eigen::VectorXd;

EnergyFunction::EnergyFunction(RadialSurface& surface,
                               std::vector<Vector3d>& extendedPrevCurve,
                               std::vector<Vector3d>& normals,
                               std::vector<Vector3d>& binormals,
                               ShellParams& parameters,
                               double radialDist,
                               std::string outputDirectory) :
                               m_surface(surface) , 
                               m_prevCurve(extendedPrevCurve) ,
                               m_normals(normals) ,
                               m_binormals(binormals) ,
                               m_parameters(parameters) ,
                               m_radialDist(radialDist),
                               m_firstRun(true),
                               m_outDirectory(outputDirectory) {
};

EnergyFunction::~EnergyFunction() = default;

double EnergyFunction::operator()(const VectorXd& inputs, VectorXd& derivatives){
    //Pre-processing

    if (m_parameters.saveEveryFrame) {
        if (m_surface.getCurveCount() > 2) {
            std::string fileName = std::to_string(m_surface.getIterCount());
            std::string outfileName = m_outDirectory+"/frame" + fileName + "surface.txt";
            std::ofstream outfile(outfileName);
            outfile << m_surface;
            outfile.close();
        }
    }
    m_surface.increaseIterCount();

    std::vector<Vector3d> nextCurve;
    int curveCount = m_surface.getCurveCount();
    int curveSize = inputs.size();
    for (int i=0; i<curveSize; i++) {
        nextCurve.push_back(m_prevCurve[i] + m_parameters.extensionLength * m_normals[i] + m_parameters.extensionLength * inputs[i] * m_binormals[i]);
        derivatives[i] = 0;
    }
    RadialSurface nextSurface(m_surface);
    nextSurface.addCurve(nextCurve);
    Vector3d zeroVec(0,0,0);

    //Use these booleans to quickly enable and disable wanted energy contributions
    
    bool lengthEnergy = false;
    bool bendEnergy = true;
    bool springEnergy = false;

    double totalEnergy = 0;
    double totalLength = 0;
    double totalBendingEnergy = 0;
    double totalSpringEnergy = 0;

    VectorXd lengthDerivs(curveSize);
    VectorXd springDerivs(curveSize);
    VectorXd bendDerivs(curveSize);

    Vector3d prev2Vec, prevVec, currentVec, nextVec, next2Vec;
    prev2Vec = nextCurve[curveSize-3];
    prevVec = nextCurve[curveSize-2];
    currentVec = nextCurve[curveSize-1];
    nextVec = nextCurve[0];
    next2Vec = nextCurve[1];
    double rescaleTerm = 1;//rescaleEnergyFunction(m_radialDist, m_parameters.radius, nextCurveLength);
    // Loop over every vertex of the new curve to find the length of the next curve
    // Also calculate the derivatives for the lengths and adds them
    
    double springNaturalLength = lengthFunction(m_radialDist,m_parameters.radius)/double(curveSize);
    for (int outerCurveIndex = 0; outerCurveIndex < curveSize; outerCurveIndex++) {
            Vector3d currentDeriv = m_parameters.extensionLength*m_binormals[outerCurveIndex];

            //Get the 5 vectors contributing to length and its derivatives, corrected for the boundaries of the curve
            prev2Vec = prevVec;
            prevVec = currentVec;
            currentVec = nextVec;
            nextVec = next2Vec;
            if (outerCurveIndex >= curveSize-2) {
                next2Vec = nextCurve[outerCurveIndex + 2 - curveSize];
            } else {
                next2Vec = nextCurve[outerCurveIndex + 2];
            }

            //Length of currently viewed edge, given by indices i and i-1
            double prevLength = (prevVec-currentVec).norm();
            //Bending contributions
            std::cout << m_parameters.bendingStiffness * bendingEnergy(prevVec, currentVec, nextVec) << std::endl;
            totalBendingEnergy += m_parameters.bendingStiffness * bendingEnergy(prevVec, currentVec, nextVec);
            bendDerivs[outerCurveIndex] += (bendEnergy) ? m_parameters.bendingStiffness * (bendingEnergyDeriv(currentVec,nextVec,next2Vec,currentDeriv,zeroVec,zeroVec) + bendingEnergyDeriv(prevVec,currentVec,nextVec,zeroVec,currentDeriv,zeroVec) + bendingEnergyDeriv(prev2Vec,prevVec,currentVec,zeroVec,zeroVec,currentDeriv)) : 0;
            
            //Length contributions
            totalLength += prevLength;
            lengthDerivs[outerCurveIndex] = (lengthEnergy) ? normDiffDeriv(prevVec, currentVec, zeroVec, currentDeriv) + normDiffDeriv(nextVec, currentVec, zeroVec, currentDeriv) : 0;
            
            //Spring Contributions
            totalSpringEnergy += 0.5 * m_parameters.strainCoeff * std::pow(prevLength - springNaturalLength,2);
            springDerivs[outerCurveIndex] = (springEnergy) ? m_parameters.strainCoeff * (normDiffDeriv(prevVec, currentVec, zeroVec, currentDeriv) * (prevLength-springNaturalLength) + normDiffDeriv(nextVec, currentVec, zeroVec, currentDeriv) * ((nextVec - currentVec).norm()-springNaturalLength)) : 0;
        }
    // Now we have the derivatives of the length at each vertex, we can multiply them by the coefficient depending on the total length to get the correct values
    double lengthDerivCoeff = m_parameters.lengthStiffness * (totalLength - lengthFunction(m_radialDist,m_parameters.radius));
    lengthDerivs *= lengthDerivCoeff;
    if (lengthEnergy) {
        derivatives += lengthDerivs;
        totalEnergy += 0.5 * m_parameters.lengthStiffness * std::pow(totalLength - lengthFunction(m_radialDist,m_parameters.radius),2);
    }
    if (springEnergy) {
        derivatives += springDerivs;
        totalEnergy += totalSpringEnergy;
    }
    if (bendEnergy) {
        derivatives += 0.5*bendDerivs;
        totalEnergy += totalBendingEnergy;
    }

    // End of loop over the previous curve vertices
    m_firstRun = false;

    derivatives *= rescaleTerm;
    return rescaleTerm * totalEnergy;
};

Vector3d EnergyFunction::normalVecDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db){
    // Returns the derivative of (a-b)/|a-b|
    double norm = (a-b).norm();
    if (norm == 0) {
        std::cout << "Division by 0!" << std::endl;
        return Vector3d(0,0,0);
    }
    return (da-db)/norm - ((da-db).dot(a-b)/std::pow(norm,3)) * (a-b);
};

double EnergyFunction::normDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db){
    // Returns the derivative of |a-b|
    double norm = (a-b).norm();
    if (norm == 0) {
        std::cout << "Division by 0!" << std::endl;
        return 0;
    }
    return (a-b).dot(da-db)/norm;
}
Vector3d EnergyFunction::normVecDeriv(Vector3d &a, Vector3d &da)
{
    // Derivative of a/|a|
    double anorm = a.norm();
    if (anorm == 0) {
        std::cout << "Division by 0!" << std::endl;
        return Vector3d(0,0,0);
    }
    double dnorm = a.dot(da)/anorm;
    return (anorm * da - dnorm * a)/std::pow(anorm,2);
}

Vector3d EnergyFunction::crossProd(Vector3d &a, Vector3d &b)
{
    return Vector3d(a[1]*b[2] - a[2]*b[1] ,a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

double EnergyFunction::dotDeriv(Vector3d &a, Vector3d &b, Vector3d &da, Vector3d &db)
{
    // Derivative of the dot product (d(a.b) = a'.b + a.b')
    return a.dot(db) + b.dot(da);
}

Vector3d EnergyFunction::crossDeriv(Vector3d &a, Vector3d &b, Vector3d &da, Vector3d &db)
{
    // Derivative of the cross product of two vectors
    return crossProd(da,b) + crossProd(a,db);
}

double EnergyFunction::lengthFunction(double t, double t0){
    double desCurv = m_parameters.desiredCurvature;
    if (desCurv < 0) {
        desCurv *= -1;
    }
    double sqrtDC = std::sqrt(desCurv);
    double length = 2 * M_PI * ( 1/sqrtDC * std::sinh(sqrtDC * (t-t0)) + t0);
    return length;
};

double EnergyFunction::bendingEnergy(Vector3d a, Vector3d b, Vector3d c){
    // 1/(l1 + l2) * tan^2 ((pi-theta)/2), theta is the angle between the c-b and a-b
    double cosAngle = (c-b).normalized().dot((a-b).normalized());
    if ((c-b).norm() == 0 && (a-b).norm() == 0) {
        std::cout << "Division by 0 error!" << std::endl;
        return 0;
    }
    if (cosAngle > 1) {
        cosAngle = 1;
    } else if (cosAngle < -1) {
        cosAngle = -1;
    }
    return (1/((c-b).norm() + (a-b).norm()) * std::pow(std::tan((M_PI - std::acos(cosAngle))/2),2));
};

double EnergyFunction::bendingEnergyDeriv(Vector3d a, Vector3d b, Vector3d c, Vector3d da, Vector3d db, Vector3d dc)
{
    // d/dx (a-b).(c-b)/|a-b||c-b|
    double dCosAngle = ((c-b).normalized().dot(normalVecDiffDeriv(a,b,da,db)) + (a-b).normalized().dot(normalVecDiffDeriv(c,b,dc,db)));

    //(a-b).(c-b)/|a-b||c-b|
    double cosAngle = (c-b).normalized().dot((a-b).normalized());

    // |a-b|
    double norm1 = (a-b).norm();

    // |c-b|
    double norm2 = (c-b).norm();

    // d/dx |a-b|
    double dnorm1 = normDiffDeriv(a,b,da,db);

    // d/dx |c-b|
    double dnorm2 = normDiffDeriv(c,b,dc,db);

    // d/dx 1/(l1+l2) tan^2((pi - std::acos(cosAngle))/2)
    double fullDeriv = ((std::pow(cosAngle,2) - 1)*(dnorm1 + dnorm2) + 2*(norm1+norm2)*dCosAngle)/(std::pow((cosAngle-1)*(norm1+norm2),2));
    return fullDeriv;
};

double EnergyFunction::rescaleEnergyFunction(double current, double init, int nextCurveLength){
    return 1/(std::pow(nextCurveLength,2));
};

// MEAN, GAUSSIAN AND AREA ENERGY CALCULATIONS, MOVED HERE TO AVOID CLUTTER IN FUNCTIONS AND NOT CURRENTLY USED

// bool areaEnergy = false;
// bool meanCurvatureEnergy = false;
// bool gaussCurvatureEnergy = false;

// // Pre-initialisations to avoid doing it every loop
//     Vector3d edge, triangleEdge1, triangleEdge2, triangleEdge3, triangleEdge4, norm1, norm2, prev2Vec, prevVec, currentVec, nextVec, next2Vec;
//     double length, edge1Length, edge2Length, edge3Length, vertexAngle;

// if (meanCurvatureEnergy || gaussCurvatureEnergy || areaEnergy) {
//         for (int curveLoc = 0; curveLoc < nextSurface.getCurveSize(curveCount-1); curveLoc++) {
//             int index = nextSurface.curveStartIndex(curveCount-1) + curveLoc;
            
//             //Code to find pairs of adjacent triangles coming off the vertex

//             auto key_selector = [](auto pair){return pair.first;};
//             double curvatureSum = 0;
//             std::vector<int> completedEdges;
//             Vector3d currentPoint = nextSurface.getPos(index);
//             std::unordered_map<int, std::vector<Triangle>> pairTriangles;
//             for (Triangle triangle : nextSurface.getNeighbourTriangles(index)) {
//                 if (triangle.vertex1 == index) {
//                     pairTriangles[triangle.vertex2].push_back(triangle);
//                     pairTriangles[triangle.vertex3].push_back(triangle);
//                 } else if (triangle.vertex2 == index) {
//                     pairTriangles[triangle.vertex1].push_back(triangle);
//                     pairTriangles[triangle.vertex3].push_back(triangle);
//                 } else {
//                     pairTriangles[triangle.vertex1].push_back(triangle);
//                     pairTriangles[triangle.vertex2].push_back(triangle);
//                 }
//             }
//             std::vector<int> keys(pairTriangles.size());
//             transform(pairTriangles.begin(), pairTriangles.end(), keys.begin(), key_selector);

//             // End of pair finding code

//             double area = 0;
//             double angleSum = 0;
            
//             // Goes through pairs of edges adjacent to eachother and calculates curvatures and areas of the triangles formed

//             for (int trianglePairIndex : keys) {
//                 if (pairTriangles[trianglePairIndex].size() == 2) {
//                     Triangle triangle1 = pairTriangles[trianglePairIndex][0];
//                     Triangle triangle2 = pairTriangles[trianglePairIndex][1];
//                     if (triangle1.vertex2 == triangle2.vertex3) {
//                         triangle1 = triangle2;
//                         triangle2 = pairTriangles[trianglePairIndex][0];
//                     }
//                     edge = nextSurface.getPos(trianglePairIndex) - currentPoint;
//                     length = edge.norm();
//                     triangleEdge1 = nextSurface.getPos(triangle1.vertex2) - currentPoint;
//                     triangleEdge2 = nextSurface.getPos(triangle1.vertex3) - currentPoint;
//                     edge1Length = (triangleEdge1).norm();
//                     edge2Length = (triangleEdge2).norm();
//                     vertexAngle = triangleEdge1.normalized().dot(triangleEdge2.normalized());
//                     triangleEdge3 = nextSurface.getPos(triangle2.vertex2) - currentPoint;
//                     triangleEdge4 = nextSurface.getPos(triangle2.vertex3) - currentPoint;
//                     norm1 = crossProd(triangleEdge1, triangleEdge2).normalized();
//                     norm2 = crossProd(triangleEdge3, triangleEdge4).normalized();
//                     area += 0.125 * edge1Length * edge2Length * std::sqrt(1-std::pow(vertexAngle,2));
//                     if (vertexAngle > 1) {
//                         vertexAngle = 1;
//                     } else if (vertexAngle < -1) {
//                         vertexAngle = -1;
//                     }
//                     angleSum += std::acos(vertexAngle); 
//                     curvatureSum += length * std::atan2(edge.normalized().dot(crossProd(norm2,norm1)), norm2.dot(norm1));
//                 }
//             }

//             if (m_firstRun) {
//                 m_origTriangleSizes.push_back(area);
//             }

//             double meanCurvature = (0.25 * curvatureSum)/area;
//             double gaussCurvature = (2 * M_PI - angleSum)/area;

//             // End of curvature and area calculations

//             // Goes through the same pairs to find the derivative of the energy function at each vertex

//             for (int trianglePairIndex : keys) {
//                 if (pairTriangles[trianglePairIndex].size() == 2) {
//                 Triangle triangle1 = pairTriangles[trianglePairIndex][0];
//                 Triangle triangle2 = pairTriangles[trianglePairIndex][1];

//                 // std::cout << "Triangle Pair" << std::endl;
//                 // std::cout << triangle1.vertex1 << " " << triangle1.vertex2 <<" " << triangle1.vertex3 << std::endl;
//                 // std::cout << triangle2.vertex1 <<" " << triangle2.vertex2 <<" " << triangle2.vertex3 << std::endl;

//                 if (triangle1.vertex2 == triangle2.vertex3) {
//                     triangle1 = triangle2;
//                     triangle2 = pairTriangles[trianglePairIndex][0];
//                 }
//                 // Find derivatives of energy contributions
//                 int nextCurveStartPoint = nextSurface.curveStartIndex(curveCount);
//                 std::vector<int> neededDerivativeIndices;
//                 if (triangle1.vertex3 >= nextCurveStartPoint) {
//                     neededDerivativeIndices.push_back(triangle1.vertex3);
//                 }
//                 if (triangle1.vertex2 >= nextCurveStartPoint) {
//                     neededDerivativeIndices.push_back(triangle1.vertex2);
//                 }
//                 if (triangle2.vertex3 >= nextCurveStartPoint) {
//                     neededDerivativeIndices.push_back(triangle2.vertex3);
//                 }

//                 for (int nextFocusPoint : neededDerivativeIndices) {
//                     // Reinitialise necessary vectors and doubles

//                     bool focusedIndex = (triangle1.vertex3 == nextFocusPoint);
//                     bool leftFocus = (triangle1.vertex2 == nextFocusPoint);
//                     int derivativeIndex = nextSurface.correctIndex(curveCount, nextFocusPoint);
                    
//                     triangleEdge1 = nextSurface.getPos(triangle1.vertex2) - currentPoint;
//                     triangleEdge2 = nextSurface.getPos(triangle1.vertex3) - currentPoint;
//                     triangleEdge3 = nextSurface.getPos(triangle2.vertex2) - currentPoint;
//                     triangleEdge4 = nextSurface.getPos(triangle2.vertex3) - currentPoint;
//                     edge1Length = (triangleEdge1).norm();
//                     edge2Length = (triangleEdge2).norm();
//                     edge3Length = (triangleEdge4).norm();
//                     norm1 = crossProd(triangleEdge1, triangleEdge2).normalized();
//                     norm2 = crossProd(triangleEdge3, triangleEdge4).normalized();

//                     // Gauss derivative calculation
                    
//                     Vector3d centrePoint = nextSurface.getPos(triangle1.vertex3);  // p2
//                     Vector3d leftNeighbour = nextSurface.getPos(triangle1.vertex2); // p1
//                     Vector3d rightNeighbour = nextSurface.getPos(triangle2.vertex3); // p3
//                     Vector3d dCentrePoint = (focusedIndex) ? m_parameters.extensionLength * m_binormals[derivativeIndex] : zeroVec; //dp2
//                     Vector3d dLeftNeighbour = (!focusedIndex && leftFocus) ? m_parameters.extensionLength * m_binormals[derivativeIndex] : zeroVec;  //dp1
//                     Vector3d dRightNeighbour = (!focusedIndex && !leftFocus) ? m_parameters.extensionLength * m_binormals[derivativeIndex] : zeroVec;   //dp3

//                     //        l2
//                     //     o------o
//                     //  l1/ p2      p3
//                     //   /
//                     //  o p1       I = (p3-p2)/|p3-p2| . (p1-p2)/|p1-p2|
//                     //             theta = acos(I)
//                     //             l1 = |p1-p2|
//                     //             l2 = |p3-p2|

//                     double vertexAngle1 = triangleEdge1.normalized().dot(triangleEdge2.normalized()); // I1
//                     double vertexAngle2 = triangleEdge3.normalized().dot(triangleEdge4.normalized()); // I2
//                     //std::cout << "Angle1 " << vertexAngle1 << " Angle2 " << vertexAngle2 << std::endl;
//                     double dt1 = angleDeriv(vertexAngle1, leftNeighbour, currentPoint, centrePoint, dLeftNeighbour, zeroVec, dCentrePoint); // I1' * acos'(I1)
//                     double dt2 = angleDeriv(vertexAngle2, centrePoint, currentPoint, rightNeighbour, dCentrePoint, zeroVec, dRightNeighbour); // I2' * acos'(I2)
//                     double dl1, dl2, dl3;
//                     double dA = triangleAreaDeriv(vertexAngle1, vertexAngle2, dt1, dt2, leftNeighbour, centrePoint, rightNeighbour, currentPoint, dLeftNeighbour, dCentrePoint, dRightNeighbour, zeroVec, dl1, dl2, dl3);
                
//                     // Mean derivative calculation
//                     double dihedralAngle;
//                     double dihedralDeriv;
//                     dihedralAngleDeriv(leftNeighbour, centrePoint, rightNeighbour, currentPoint, dLeftNeighbour, dCentrePoint, dRightNeighbour, zeroVec, dihedralAngle, dihedralDeriv);
                    
//                     double dGauss = (focusedIndex) ? -(dt1+dt2)/area - dA * (gaussCurvature)/area : 0;
//                     double dMean = (focusedIndex) ? 0.25 * (dl2 * dihedralAngle + edge2Length * dihedralDeriv)/area - dA*meanCurvature/area : 0.25 * edge2Length * dihedralDeriv/area;
//                     double dStretch = (focusedIndex) ? m_parameters.strainCoeff * (area - m_origTriangleSizes[curveLoc]) * dA : 0;
                    
//                     double meanCurvatureDerivContribution = (meanCurvatureEnergy) ? m_parameters.meanStiffness * (meanCurvature - m_parameters.desiredCurvature) * dMean : 0;
//                     double gaussCurvatureDerivContribution = (gaussCurvatureEnergy) ? - m_parameters.gaussStiffness * dGauss : 0;
//                     double stretchDerivContribution =  (areaEnergy) ? m_parameters.strainCoeff * (dStretch) : 0;

//                     // Finally add the contributions to the derivatives
//                     derivatives[derivativeIndex] +=  meanCurvatureDerivContribution + gaussCurvatureDerivContribution + stretchDerivContribution;
//                     }
//                 }
//             }
//             double meanCurvatureContribution = (meanCurvatureEnergy) ? m_parameters.meanStiffness * 0.5 * std::pow(meanCurvature - m_parameters.desiredCurvature,2) : 0;
//             double gaussCurvatureContribution = (gaussCurvatureEnergy) ? - m_parameters.gaussStiffness * gaussCurvature : 0;
//             double stretchingEnergyContribution = (areaEnergy) ? m_parameters.strainCoeff * 0.5 * std::pow(area - m_origTriangleSizes[curveLoc],2) : 0;

//             // Add on the contribution to enery from each vertex
//             totalEnergy += meanCurvatureContribution + gaussCurvatureContribution + stretchingEnergyContribution;
//         }
//     }


// UNUSED FUNCTIONS STORED DOWN HERE

// DIHEDRAL ANGLE DERIVATIVE FUNCTION

void EnergyFunction::dihedralAngleDeriv(Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &d, Vector3d &da, Vector3d &db, Vector3d &dc, Vector3d &dd, double& angleResult, double& derivResult)
{
    // Finds the dihedral angle derivative across the edge b-d, where d is the centre vertex and a,b,c are the points that make up the faces
    Vector3d edge1 = a-d;
    Vector3d edge2 = b-d;
    Vector3d edge3 = c-d;
    Vector3d dedge1 = da-dd;
    Vector3d dedge2 = db-dd;
    Vector3d dedge3 = dc-dd;
    Vector3d edge2Norm = edge2.normalized();
    Vector3d dedge2Norm = normVecDeriv(edge2, dedge2);
    Vector3d unNormed1 = crossProd(edge1, edge2);
    Vector3d unNormed2 = crossProd(edge2, edge3);
    Vector3d dUnNormed1 = crossDeriv(edge1, edge2, dedge1, dedge2);
    Vector3d dUnNormed2 = crossDeriv(edge2, edge3, dedge2, dedge3);
    Vector3d norm1 = unNormed1.normalized();
    Vector3d norm2 = unNormed2.normalized();
    Vector3d dNorm1 = normVecDeriv(unNormed1, dUnNormed1);
    Vector3d dNorm2 = normVecDeriv(unNormed2, dUnNormed2);
    Vector3d crossNorm = crossProd(norm1, norm2);
    Vector3d dCrossNorm = crossDeriv(norm1, norm2, dNorm1, dNorm2);
    double atan1 = edge2Norm.dot(crossProd(norm1, norm2));
    double atan2 = norm2.dot(norm1);
    double datan1 = dotDeriv(edge2Norm, crossNorm, dedge2Norm, dCrossNorm);
    double datan2 = dotDeriv(norm2, norm1, dNorm2, dNorm1);
    angleResult = std::atan2(atan1, atan2);
    derivResult = atan2Deriv(atan1, atan2, datan1, datan2);
}

double EnergyFunction::triangleAreaDeriv(double insideAngle1, double insideAngle2, double dangle1, double dangle2 ,Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &d, Vector3d &da, Vector3d &db, Vector3d &dc, Vector3d &dd, double& dl1, double& dl2, double& dl3)
{
    // Takes the point and derivatives of the only contributing terms to the area of a triangle face
    // CurrentPoint is d, others are a,b,c
    double l1 = (a-d).norm();
    double l2 = (b-d).norm();
    double l3 = (c-d).norm();
    dl1 = normDiffDeriv(a,d,da,dd);
    dl2 = normDiffDeriv(b,d,db,dd);
    dl3 = normDiffDeriv(c,d,dc,dd);
    double part1 = 0.125 * (l1 * l2 * dangle1 * insideAngle1 + (l1*dl2 + l2*dl1)*std::sqrt(1-std::pow(insideAngle1,2)));
    double part2 = 0.125 * (l2 * l3 * dangle2 * insideAngle2 + (l2*dl3 + l3*dl2)*std::sqrt(1-std::pow(insideAngle2,2)));
    return part1 + part2;
}

double EnergyFunction::angleDeriv(double insideAngle, Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &da, Vector3d &db, Vector3d &dc)
{
    // Finds the derivative for the angle between the vectors (a-b) and (c-b)
    
    Vector3d normalVecDeriv1 = normalVecDiffDeriv(a,b,da,db);
    Vector3d normalVecDeriv2 = normalVecDiffDeriv(c,b,dc,db);
    Vector3d normalisedVec1 = (a-b).normalized();
    Vector3d normalisedVec2 = (c-b).normalized();
    double insideAngleDeriv = dotDeriv(normalisedVec1, normalisedVec2, normalVecDeriv1, normalVecDeriv2);
    if (insideAngleDeriv > 1) {
        std::cout << "Division by 0!" << std::endl;
        insideAngleDeriv = 1;
        return 0;
    } else if (insideAngleDeriv < -1) {
        std::cout << "Division by 0!" << std::endl;
        insideAngleDeriv = -1;
        return 0;
    }

    return insideAngleDeriv * (-1/std::sqrt(1-std::pow(insideAngle,2)));
}

double EnergyFunction::atan2Deriv(double x, double y, double dx, double dy)
{
    // Returns the derivative of atan2(x,y)
    return (x*dy - y*dx)/(std::pow(x,2) + std::pow(y,2));
};