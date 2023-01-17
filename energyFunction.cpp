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
                               m_radialDist(radialDist) {
};

EnergyFunction::~EnergyFunction() = default;

double EnergyFunction::operator()(const VectorXd& inputs, VectorXd& derivatives){
    std::vector<Vector3d> nextCurve;
    int curveCount = m_surface.getCurveCount();
    double sParam;
    int curveSize = inputs.size();
    for (int i=0; i<curveSize; i++) {
        nextCurve.push_back(m_prevCurve[i] + m_parameters.extensionLength * m_normals[i] + m_parameters.extensionLength * inputs[i] * m_binormals[i]);
        derivatives[i] = 0;
    }
    RadialSurface nextSurface(m_surface);
    nextSurface.addCurve(nextCurve);
    
    //Reimplementation of finding curvature and derivatives, moved into this space to keep it all together

    // Pre-initialisations to avoid doing it every loop
    Vector3d edge;
    double length;
    Vector3d triangleEdge1;
    Vector3d triangleEdge2;
    double edge1Length;
    double edge2Length;
    double edge3Length;
    double vertexAngle;
    Vector3d triangleEdge3;
    Vector3d triangleEdge4;
    Vector3d norm1;
    Vector3d norm2;
    // End
    double totalEnergy = 0;
    for (int curveLoc = 0; curveLoc < nextSurface.getCurveLength(curveCount-1); curveLoc++) {
        int index = nextSurface.curveStartIndex(curveCount-1) + curveLoc;
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
        double area = 0;
        double angleSum = 0;
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
        double meanCurvature = (0.25 * curvatureSum)/area;
        double gaussCurvature = (2 * M_PI - angleSum)/area;

        for (int trianglePairIndex : keys) {
            if (pairTriangles[trianglePairIndex].size() == 2) {
            Triangle triangle1 = pairTriangles[trianglePairIndex][0];
            Triangle triangle2 = pairTriangles[trianglePairIndex][1];
            if (triangle1.vertex2 == triangle2.vertex3) {
                triangle1 = triangle2;
                triangle2 = pairTriangles[trianglePairIndex][0];
            }
            // Find derivatives of energy contributions
            if (triangle1.vertex3 >= nextSurface.curveStartIndex(curveCount)) {

                // Reinitialise necessary vectors and doubles
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
                int derivativeSurfaceIndex = triangle1.vertex3;
                int derivativeIndex = nextSurface.correctIndex(curveCount, triangle1.vertex3);
                Vector3d derivPoint = nextSurface.getPos(derivativeSurfaceIndex);
                Vector3d derivPointNeighbour1 = nextSurface.getPos(triangle1.vertex2);
                Vector3d derivPointNeighbour2 = nextSurface.getPos(triangle2.vertex3);
                Vector3d dDerivPoint = m_parameters.extensionLength * m_binormals[derivativeIndex];
                Vector3d zeroVec(0,0,0);
                double vertexAngle1 = triangleEdge1.normalized().dot(triangleEdge2.normalized());
                double vertexAngle2 = triangleEdge3.normalized().dot(triangleEdge4.normalized());
                double dt1 = angleDeriv(vertexAngle1, derivPointNeighbour1, currentPoint, derivPoint, zeroVec, zeroVec, dDerivPoint);
                double dt2 = angleDeriv(vertexAngle2, derivPoint, currentPoint, derivPointNeighbour2, dDerivPoint, zeroVec, zeroVec);
                double dA = triangleAreaDeriv(vertexAngle1, vertexAngle2, dt1, dt2, derivPointNeighbour1, derivPoint, derivPointNeighbour2, currentPoint, zeroVec, dDerivPoint, zeroVec, zeroVec);
                double dGauss = -(dt1+dt2)/area - dA * (gaussCurvature)/std::pow(area,2);
                
                // Mean derivative calculation
                Vector3d normalisedEdge1 = triangleEdge1.normalized();
                Vector3d normalisedEdge2 = triangleEdge3.normalized();
                Vector3d normalisedEdge3 = triangleEdge4.normalized();
                Vector3d edgeVecDeriv = normVecDeriv(triangleEdge2, dDerivPoint);
                Vector3d norm1Numerator = crossProd(triangleEdge1, triangleEdge2);
                Vector3d norm2Numerator = crossProd(triangleEdge3, triangleEdge4);
                Vector3d dnorm1Numerator = crossDeriv(triangleEdge1, triangleEdge2, zeroVec, edgeVecDeriv);
                Vector3d dnorm2Numerator = crossDeriv(triangleEdge3, triangleEdge4, edgeVecDeriv, zeroVec);
                Vector3d dnorm1 = normVecDeriv(norm1Numerator, dnorm1Numerator);
                Vector3d dnorm2 = normVecDeriv(norm2Numerator, dnorm2Numerator);
                Vector3d crossNorms = crossProd(norm1, norm2);
                double atanp1 = normalisedEdge2.dot(crossNorms);
                double atanp2 = norm2.dot(norm1);
                Vector3d insideCrossDeriv = crossDeriv(norm1, norm2, dnorm1, dnorm2);
                double datanp1 = dotDeriv(normalisedEdge2, crossNorms, edgeVecDeriv, insideCrossDeriv);
                double datanp2 = dotDeriv(norm2, norm1, dnorm2, dnorm1);
                double dda = atan2Deriv(atanp1, atanp2, datanp1, datanp2);
                double dl1 = normDiffDeriv(derivPoint, currentPoint, dDerivPoint, zeroVec);
                double dihedralAngle = std::atan2(edge.normalized().dot(crossProd(norm2,norm1)), norm2.dot(norm1));
                double dMean = 0.25 * ((dl1 * dihedralAngle + edge2Length * dda)/area - dA*edge2Length*dihedralAngle/std::pow(area,2));

                // Each point has 3 contributions, so we calculate those contributions on the neighbouring vertices and add them
                if (triangle1.vertex2 >= nextSurface.curveStartIndex(curveCount)){
                    int derivativeIndexNghb1 = nextSurface.correctIndex(curveCount, triangle1.vertex2);
                    Vector3d dDerivPointNeighbour1 =  m_parameters.extensionLength * m_binormals[derivativeIndexNghb1];
                    Vector3d nbghEdgeVecDeriv1 = normalVecDiffDeriv(derivPointNeighbour1, currentPoint, dDerivPointNeighbour1, zeroVec);
                    Vector3d dnbghNorm1Numerator = crossDeriv(triangleEdge1, triangleEdge2, nbghEdgeVecDeriv1, zeroVec);
                    Vector3d dNghbNorm1 = normVecDeriv(norm1Numerator, dnbghNorm1Numerator);
                    Vector3d nghbInsideCrossDeriv1 = crossDeriv(norm1, norm2, dNghbNorm1, zeroVec);
                    double dNghbAtanp11 = dotDeriv(normalisedEdge2, crossNorms, zeroVec, nghbInsideCrossDeriv1);
                    double dNghbAtanp21 = dotDeriv(norm2, norm1, zeroVec, dNghbNorm1);
                    double ddNghbA1 = atan2Deriv(atanp1, atanp2, dNghbAtanp11, dNghbAtanp21);
                    double nghbDt11 = angleDeriv(vertexAngle1, derivPointNeighbour1, currentPoint, derivPoint, dDerivPointNeighbour1, zeroVec, zeroVec);
                    double nghbDt21 = angleDeriv(vertexAngle2, derivPoint, currentPoint, derivPointNeighbour2, zeroVec, zeroVec, zeroVec);
                    double nghbDA1 = triangleAreaDeriv(vertexAngle1, vertexAngle2, nghbDt11, nghbDt21, derivPointNeighbour1, derivPoint, derivPointNeighbour2, currentPoint, dDerivPointNeighbour1, zeroVec, zeroVec, zeroVec);
                    derivatives[derivativeIndexNghb1] += m_parameters.stiffnessRatio * (meanCurvature - m_parameters.desiredCurvature) * 0.25 * edge1Length * (ddNghbA1/area - dihedralAngle*nghbDA1/std::pow(area,2));
                    
                }

                if (triangle2.vertex3 >= nextSurface.curveStartIndex(curveCount)){
                    int derivativeIndexNghb2 = nextSurface.correctIndex(curveCount, triangle2.vertex3);
                    Vector3d dDerivPointNeighbour2 =  m_parameters.extensionLength * m_binormals[derivativeIndexNghb2];
                    Vector3d nbghEdgeVecDeriv2 = normalVecDiffDeriv(derivPointNeighbour2, currentPoint, dDerivPointNeighbour2, zeroVec);
                    Vector3d dnbghNorm2Numerator = crossDeriv(triangleEdge3, triangleEdge4, zeroVec, nbghEdgeVecDeriv2);
                    Vector3d dNghbNorm2 = normVecDeriv(norm2Numerator, dnbghNorm2Numerator);
                    Vector3d nghbInsideCrossDeriv2 = crossDeriv(norm1, norm2, zeroVec, dNghbNorm2);
                    double dNghbAtanp12 = dotDeriv(normalisedEdge2, crossNorms, zeroVec, nghbInsideCrossDeriv2);
                    double dNghbAtanp22 = dotDeriv(norm2, norm1, dNghbNorm2, zeroVec);
                    double ddNghbA2 = atan2Deriv(atanp1, atanp2, dNghbAtanp12, dNghbAtanp22);
                    double nghbDt12 = angleDeriv(vertexAngle1, derivPointNeighbour1, currentPoint, derivPoint, zeroVec, zeroVec, zeroVec);
                    double nghbDt22 = angleDeriv(vertexAngle2, derivPoint, currentPoint, derivPointNeighbour2, zeroVec, zeroVec, dDerivPointNeighbour2);
                    double nghbDA2 = triangleAreaDeriv(vertexAngle1, vertexAngle2, nghbDt12, nghbDt22, derivPointNeighbour1, derivPoint, derivPointNeighbour2, currentPoint, zeroVec, zeroVec, dDerivPointNeighbour2, zeroVec);
                    derivatives[derivativeIndexNghb2] += m_parameters.stiffnessRatio * (meanCurvature - m_parameters.desiredCurvature) * 0.25 * edge3Length * (ddNghbA2/area - dihedralAngle*nghbDA2/std::pow(area,2));
                    
                }
                // Finally add the contributions to the derivatives
                
                derivatives[derivativeIndex] += m_parameters.stiffnessRatio * (meanCurvature - m_parameters.desiredCurvature) * dMean + (gaussCurvature - m_parameters.desiredCurvature)* dGauss; 
                
                }
            }
        }
        totalEnergy += m_parameters.stiffnessRatio * 0.5 * std::pow(meanCurvature - m_parameters.desiredCurvature,2) + 0.5 * std::pow(gaussCurvature - m_parameters.desiredCurvature,2);
    }
    return totalEnergy;
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
        totalMean += m_parameters.stiffnessRatio / 2 * std::pow(currentMean - m_parameters.desiredCurvature,2);
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

double EnergyFunction::triangleAreaDeriv(double angle1, double angle2, double dangle1, double dangle2 ,Vector3d &a, Vector3d &b, Vector3d &c, Vector3d &d, Vector3d &da, Vector3d &db, Vector3d &dc, Vector3d &dd)
{
    // Takes the point and derivatives of the only contributing terms to the area of a triangle face
    // CurrentPoint is d, others are a,b,c
    double l1 = (a-d).norm();
    double l2 = (b-d).norm();
    double l3 = (c-d).norm();
    double dl1 = normDiffDeriv(a,d,da,dd);
    double dl2 = normDiffDeriv(b,d,db,dd);
    double dl3 = normDiffDeriv(c,d,dc,dd);
    double part1 = 0.125 * (l1 * l2 * dangle1 * std::cos(angle1) + (l1*dl2 + l2*dl1)*std::sin(angle1));
    double part2 = 0.125 * (l2 * l3 * dangle2 * std::cos(angle2) + (l2*dl3 + l3*dl2)*std::sin(angle2));
    return part1 + part2;
}
