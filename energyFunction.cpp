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
                               double radialDist) :
                               m_surface(surface) , 
                               m_prevCurve(extendedPrevCurve) ,
                               m_normals(normals) ,
                               m_binormals(binormals) ,
                               m_parameters(parameters) ,
                               m_radialDist(radialDist),
                               m_firstRun(true) {
};

EnergyFunction::~EnergyFunction() = default;

double EnergyFunction::operator()(const VectorXd& inputs, VectorXd& derivatives){

    
    // Current energy function is defined by 
    //
    //  B/2 * (H-c0)**2 + K + C/2 * (A-A0)^2 + SL/2 * (L-L0) ** 2
    //
    // Summed over all vertices of the previous ring
    //
    // H: Mean curvature, K: Gauss curvature, B: Bending stiffness, C: Stretch Resistance, A: Current triangle area, A0: Original vertex area

    //Use these booleans to quickly enable and disable wanted energy contributions
    bool meanCurvatureEnergy = true;
    bool gaussCurvatureEnergy = true;
    bool lengthEnergy = false;
    bool areaEnergy = false;

    std::vector<Vector3d> nextCurve;
    int curveCount = m_surface.getCurveCount();
    int curveSize = inputs.size();
    for (int i=0; i<curveSize; i++) {
        nextCurve.push_back(m_prevCurve[i] + m_parameters.extensionLength * m_normals[i] + m_parameters.extensionLength * inputs[i] * m_binormals[i]);
        derivatives[i] = 0;
    }
    RadialSurface nextSurface(m_surface);
    nextSurface.addCurve(nextCurve);
    // Pre-initialisations to avoid doing it every loop
    Vector3d edge, triangleEdge1, triangleEdge2, triangleEdge3, triangleEdge4, norm1, norm2, prevVec, currentVec, nextVec;
    double length, edge1Length, edge2Length, edge3Length, vertexAngle;
    int nextCurveLength = nextSurface.getCurveLength(curveCount);
    double totalLength = 0;

    double totalEnergy = 0;
    double rescaleTerm = rescaleEnergyFunction(m_radialDist, m_parameters.radius, nextCurveLength);
    // Loop over every vertex of the new curve to find the length of the next curve
    // Also calculate the derivatives for the lengths and adds them
    
    for (int outerCurveIndex = 0; outerCurveIndex < nextCurveLength; outerCurveIndex++) {
        //Get the 3 vectors contributing to length and its derivatives, corrected for the boundaries of the curve
        currentVec = nextCurve[outerCurveIndex];
        if (outerCurveIndex == nextCurveLength-1) {
            nextVec = nextCurve[0];
            prevVec = nextCurve[outerCurveIndex-1];
        } else if (outerCurveIndex == 0) {
            nextVec = nextCurve[outerCurveIndex+1];
            prevVec = nextCurve[nextCurveLength-1];
        } else {
            nextVec = nextCurve[outerCurveIndex+1];
            prevVec = nextCurve[outerCurveIndex-1];
        }
        // With the vectors, calculate the length and add its contribution (only between two of the vectors as it sums over all vectors)
        totalLength += (currentVec - prevVec).norm();
        Vector3d zeroVec(0,0,0);
        Vector3d currentDeriv = m_parameters.extensionLength*m_binormals[outerCurveIndex];
        derivatives[outerCurveIndex] = (lengthEnergy) ? normDiffDeriv(prevVec, currentVec, zeroVec, currentDeriv) + normDiffDeriv(nextVec, currentVec, zeroVec, currentDeriv) : 0;
        }
    // Now we have the derivatives of the length at each vertex, we can multiply them by the coefficient depending on the total length to get the correct values
    double lengthDerivCoeff = m_parameters.lengthStiffness * (totalLength - lengthFunction(m_radialDist,m_parameters.radius));
    derivatives *= lengthDerivCoeff;

    for (int curveLoc = 0; curveLoc < nextSurface.getCurveLength(curveCount-1); curveLoc++) {
        int index = nextSurface.curveStartIndex(curveCount-1) + curveLoc;
        
        //Code to find pairs of adjacent triangles coming off the vertex

        auto key_selector = [](auto pair){return pair.first;};
        double curvatureSum = 0;
        std::vector<int> completedEdges;
        Vector3d currentPoint = nextSurface.getPos(index);
        std::unordered_map<int, std::vector<Triangle>> pairTriangles;
        for (Triangle triangle : nextSurface.getNeighbourTriangles(index)) {
            if (triangle.vertex1 == index) {
                pairTriangles[triangle.vertex2].push_back(triangle);
                pairTriangles[triangle.vertex3].push_back(triangle);
            } else if (triangle.vertex2 == index) {
                pairTriangles[triangle.vertex1].push_back(triangle);
                pairTriangles[triangle.vertex3].push_back(triangle);
            } else {
                pairTriangles[triangle.vertex1].push_back(triangle);
                pairTriangles[triangle.vertex2].push_back(triangle);
            }
        }
        std::vector<int> keys(pairTriangles.size());
        transform(pairTriangles.begin(), pairTriangles.end(), keys.begin(), key_selector);

        // End of pair finding code

        double area = 0;
        double angleSum = 0;
        
        // Goes through pairs of edges adjacent to eachother and calculates curvatures and areas of the triangles formed

        for (int trianglePairIndex : keys) {
            if (pairTriangles[trianglePairIndex].size() == 2) {
                Triangle triangle1 = pairTriangles[trianglePairIndex][0];
                Triangle triangle2 = pairTriangles[trianglePairIndex][1];
                if (triangle1.vertex2 == triangle2.vertex3) {
                    triangle1 = triangle2;
                    triangle2 = pairTriangles[trianglePairIndex][0];
                }
                edge = nextSurface.getPos(trianglePairIndex) - currentPoint;
                length = edge.norm();
                triangleEdge1 = nextSurface.getPos(triangle1.vertex2) - currentPoint;
                triangleEdge2 = nextSurface.getPos(triangle1.vertex3) - currentPoint;
                edge1Length = (triangleEdge1).norm();
                edge2Length = (triangleEdge2).norm();
                vertexAngle = std::acos(triangleEdge1.normalized().dot(triangleEdge2.normalized()));
                triangleEdge3 = nextSurface.getPos(triangle2.vertex2) - currentPoint;
                triangleEdge4 = nextSurface.getPos(triangle2.vertex3) - currentPoint;
                norm1 = crossProd(triangleEdge1, triangleEdge2).normalized();
                norm2 = crossProd(triangleEdge3, triangleEdge4).normalized();
                area += 0.125 * edge1Length * edge2Length * std::sin(vertexAngle);
                angleSum += vertexAngle; 
                curvatureSum += length * std::atan2(edge.normalized().dot(crossProd(norm2,norm1)), norm2.dot(norm1));
            }
        }

        if (m_firstRun) {
            m_origTriangleSizes.push_back(area);
        }

        double meanCurvature = (0.25 * curvatureSum)/area;
        double gaussCurvature = (2 * M_PI - angleSum)/area;

        // End of curvature and area calculations

        // Goes through the same pairs to find the derivative of the energy function at each vertex

        for (int trianglePairIndex : keys) {
            if (pairTriangles[trianglePairIndex].size() == 2) {
            Triangle triangle1 = pairTriangles[trianglePairIndex][0];
            Triangle triangle2 = pairTriangles[trianglePairIndex][1];

            // std::cout << "Triangle Pair" << std::endl;
            // std::cout << triangle1.vertex1 << " " << triangle1.vertex2 <<" " << triangle1.vertex3 << std::endl;
            // std::cout << triangle2.vertex1 <<" " << triangle2.vertex2 <<" " << triangle2.vertex3 << std::endl;

            if (triangle1.vertex2 == triangle2.vertex3) {
                triangle1 = triangle2;
                triangle2 = pairTriangles[trianglePairIndex][0];
            }
            // Find derivatives of energy contributions
            int nextCurveStartPoint = nextSurface.curveStartIndex(curveCount);
            std::vector<int> neededDerivativeIndices;
            if (triangle1.vertex3 >= nextCurveStartPoint) {
                neededDerivativeIndices.push_back(triangle1.vertex3);
            }
            if (triangle1.vertex2 >= nextCurveStartPoint) {
                neededDerivativeIndices.push_back(triangle1.vertex2);
            }
            if (triangle2.vertex3 >= nextCurveStartPoint) {
                neededDerivativeIndices.push_back(triangle2.vertex3);
            }

            for (int nextFocusPoint : neededDerivativeIndices) {
                // Reinitialise necessary vectors and doubles

                bool focusedIndex = (triangle1.vertex3 == nextFocusPoint);
                bool leftFocus = (triangle1.vertex2 == nextFocusPoint);
                int derivativeIndex = nextSurface.correctIndex(curveCount, nextFocusPoint);
                
                triangleEdge1 = nextSurface.getPos(triangle1.vertex2) - currentPoint;
                triangleEdge2 = nextSurface.getPos(triangle1.vertex3) - currentPoint;
                triangleEdge3 = nextSurface.getPos(triangle2.vertex2) - currentPoint;
                triangleEdge4 = nextSurface.getPos(triangle2.vertex3) - currentPoint;
                edge1Length = (triangleEdge1).norm();
                edge2Length = (triangleEdge2).norm();
                edge3Length = (triangleEdge4).norm();
                norm1 = crossProd(triangleEdge1, triangleEdge2).normalized();
                norm2 = crossProd(triangleEdge3, triangleEdge4).normalized();

                // Gauss derivative calculation
                
                Vector3d zeroVec(0,0,0);
                Vector3d centrePoint = nextSurface.getPos(triangle1.vertex3);
                Vector3d leftNeighbour = nextSurface.getPos(triangle1.vertex2);
                Vector3d rightNeighbour = nextSurface.getPos(triangle2.vertex3);
                Vector3d dCentrePoint = (focusedIndex) ? m_parameters.extensionLength * m_binormals[derivativeIndex] : zeroVec;
                Vector3d dLeftNeighbour = (!focusedIndex && leftFocus) ? m_parameters.extensionLength * m_binormals[derivativeIndex] : zeroVec;
                Vector3d dRightNeighbour = (!focusedIndex && !leftFocus) ? m_parameters.extensionLength * m_binormals[derivativeIndex] : zeroVec;

                double vertexAngle1 = triangleEdge1.normalized().dot(triangleEdge2.normalized());
                double vertexAngle2 = triangleEdge3.normalized().dot(triangleEdge4.normalized());
                double dt1 = angleDeriv(vertexAngle1, leftNeighbour, currentPoint, centrePoint, dLeftNeighbour, zeroVec, dCentrePoint);
                double dt2 = angleDeriv(vertexAngle2, centrePoint, currentPoint, rightNeighbour, dCentrePoint, zeroVec, dRightNeighbour);
                double dl1, dl2, dl3;
                double dA = triangleAreaDeriv(vertexAngle1, vertexAngle2, dt1, dt2, leftNeighbour, centrePoint, rightNeighbour, currentPoint, dLeftNeighbour, dCentrePoint, dRightNeighbour, zeroVec, dl1, dl2, dl3);
            
                // Mean derivative calculation
                double dihedralAngle;
                double dihedralDeriv;
                dihedralAngleDeriv(leftNeighbour, centrePoint, rightNeighbour, currentPoint, dLeftNeighbour, dCentrePoint, dRightNeighbour, zeroVec, dihedralAngle, dihedralDeriv);
                
                double dGauss = (focusedIndex) ? -(dt1+dt2)/area - dA * (gaussCurvature)/std::pow(area,2) : 0;
                // if (derivativeIndex == 10) {
                    
                //     if (focusedIndex) {
                //         std::cout << "Curve Index" << curveLoc << " " <<  dGauss << std::endl;
                //         std::cout << "Triangle Contributer" << std::endl;
                //         std::cout << triangle1.vertex1 << " " << triangle1.vertex2 <<" " << triangle1.vertex3 << std::endl;
                //         std::cout << triangle2.vertex1 <<" " << triangle2.vertex2 <<" " << triangle2.vertex3 << std::endl;
                //     }
                // } 
                double dMean = (focusedIndex) ? 0.25 * (dl2 * dihedralAngle + edge2Length * dihedralDeriv)/area - dA*meanCurvature/area : 0.25 * edge2Length * dihedralDeriv/area;
                double dStretch = (focusedIndex) ? m_parameters.strainCoeff * (area - m_origTriangleSizes[curveLoc]) * dA : 0;
                
                double meanCurvatureDerivContribution = (meanCurvatureEnergy) ? m_parameters.meanStiffness * (meanCurvature) * dMean : 0;
                double gaussCurvatureDerivContribution = (gaussCurvatureEnergy) ? m_parameters.gaussStiffness * (gaussCurvature - m_parameters.desiredCurvature) * dGauss : 0;
                double stretchDerivContribution =  (areaEnergy) ? m_parameters.strainCoeff * (dStretch) : 0;

                // Finally add the contributions to the derivatives
                derivatives[derivativeIndex] +=  meanCurvatureDerivContribution + gaussCurvatureDerivContribution + stretchDerivContribution;
                }
            }
        }
        double meanCurvatureContribution = (meanCurvatureEnergy) ? m_parameters.meanStiffness * 0.5 * std::pow(meanCurvature,2) : 0;
        double gaussCurvatureContribution = (gaussCurvatureEnergy) ? m_parameters.gaussStiffness * std::pow(gaussCurvature - m_parameters.desiredCurvature,2) : 0;
        double stretchingEnergyContribution = (areaEnergy) ? m_parameters.strainCoeff * 0.5 * std::pow(area - m_origTriangleSizes[curveLoc],2) : 0;
        

        // Add on the contribution to enery from each vertex
        totalEnergy += meanCurvatureContribution + gaussCurvatureContribution + stretchingEnergyContribution;
    }
    
    double lengthEnergyContribution = (lengthEnergy) ? 0.5 * m_parameters.lengthStiffness * std::pow(totalLength - lengthFunction(m_radialDist,m_parameters.radius),2) : 0;
    totalEnergy += lengthEnergyContribution;
    // End of loop over the previous curve vertices
    m_firstRun = false;

    derivatives *= rescaleTerm;
    return rescaleTerm * totalEnergy;
};

double EnergyFunction::evalEnergy(RadialSurface& extendedSurface)
{
    int curveCount = extendedSurface.getCurveCount();
    // Here find the mean and gaussian curvatures of the surface at every vertex
    double totalGauss = 0;
    double totalMean = 0;
    double currentGauss;
    double currentMean;
    for (int i = 0; i < extendedSurface.getCurveLength(curveCount-2); i++) {
        int index = extendedSurface.curveStartIndex(curveCount-2) + i;
        extendedSurface.curvatures(index, currentGauss, currentMean);
        totalGauss += currentGauss;
        totalMean += m_parameters.meanStiffness / 2 * std::pow(currentMean - m_parameters.desiredCurvature,2);
    }
    return totalMean + totalGauss;
};

Vector3d EnergyFunction::normalVecDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db){
    // Returns the derivative of (a-b)/|a-b|
    double norm = (a-b).norm();
    return (da-db)/norm - ((da-db).dot(a-b)/std::pow(norm,3)) * (a-b);
};

double EnergyFunction::normDiffDeriv(Vector3d& a, Vector3d& b, Vector3d& da, Vector3d& db){
    // Returns the derivative of |a-b|
    double norm = (a-b).norm();
    if (norm == 0) {
        return 0;
    }
    return (a-b).dot(da-db)/norm;
}
Vector3d EnergyFunction::normVecDeriv(Vector3d &a, Vector3d &da)
{
    // Derivative of a/|a|
    double anorm = a.norm();
    double dnorm = a.dot(da)/anorm;
    return (anorm * da - dnorm * a)/std::pow(anorm,2);
}
double EnergyFunction::atan2Deriv(double x, double y, double dx, double dy)
{
    // Returns the derivative of atan2(x,y)
    return (x*dy - y*dx)/(std::pow(x,2) + std::pow(y,2));
};

Vector3d EnergyFunction::crossProd(Vector3d &a, Vector3d &b)
{
    return Vector3d(a[1]*b[2] - a[2]*b[1] ,a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

double EnergyFunction::angleDeriv(double angle, Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &da, Vector3d &db, Vector3d &dc)
{
    // Finds the derivative for the angle between the vectors (a-b) and (c-b)
    
    Vector3d normalVecDeriv1 = normalVecDiffDeriv(a,b,da,db);
    Vector3d normalVecDeriv2 = normalVecDiffDeriv(c,b,dc,db);
    Vector3d normalisedVec1 = (a-b).normalized();
    Vector3d normalisedVec2 = (c-b).normalized();
    double insideAngleDeriv = dotDeriv(normalisedVec1, normalisedVec2, normalVecDeriv1, normalVecDeriv2);
    return insideAngleDeriv * (-1/std::sqrt(1-std::pow(angle,2)));
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

double EnergyFunction::triangleAreaDeriv(double angle1, double angle2, double dangle1, double dangle2 ,Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &d, Vector3d &da, Vector3d &db, Vector3d &dc, Vector3d &dd, double& dl1, double& dl2, double& dl3)
{
    // Takes the point and derivatives of the only contributing terms to the area of a triangle face
    // CurrentPoint is d, others are a,b,c
    double l1 = (a-d).norm();
    double l2 = (b-d).norm();
    double l3 = (c-d).norm();
    dl1 = normDiffDeriv(a,d,da,dd);
    dl2 = normDiffDeriv(b,d,db,dd);
    dl3 = normDiffDeriv(c,d,dc,dd);
    double part1 = 0.125 * (l1 * l2 * dangle1 * std::cos(angle1) + (l1*dl2 + l2*dl1)*std::sin(angle1));
    double part2 = 0.125 * (l2 * l3 * dangle2 * std::cos(angle2) + (l2*dl3 + l3*dl2)*std::sin(angle2));
    return part1 + part2;
}

double EnergyFunction::lengthFunction(double t, double t0){
    double desCurv = m_parameters.desiredCurvature;
    if (desCurv < 0) {
        desCurv *= -1;
    }
    double sqrtDC = std::sqrt(desCurv);
    return 2 * M_PI * ( 1/sqrtDC * std::sinh(sqrtDC * (t-t0)) + t0);
};

double EnergyFunction::rescaleEnergyFunction(double t, double t0, int nextCurveLength){
    //return 1/(std::pow(nextCurveLength,2) * lengthFunction(t,t0));
    double multiplier = m_parameters.lengthStiffness;
    if (m_parameters.lengthStiffness < 1) {
        multiplier = 1/m_parameters.lengthStiffness;
    }
    if (m_parameters.meanStiffness < 1) {
        multiplier *= 1/m_parameters.meanStiffness;
    } else {
        multiplier *= m_parameters.meanStiffness;
    }
    return 1/(lengthFunction(t,t0));
};