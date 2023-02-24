#include "surface.h"

Surface::Surface() : m_pointCount(0)
{
}

Surface::~Surface()
{
}

void Surface::addPoint(Vector3d& newPoint)
{   
    m_pointCount += 1;
    m_points.push_back(newPoint);
    m_faces[m_pointCount-1];
    return;
}

int Surface::surfaceSize()
{
    return m_pointCount;
}

void Surface::addTriangle(Triangle triangle)
{
    m_faces[triangle.vertex1].push_back(triangle);
    m_faces[triangle.vertex2].push_back(Triangle(triangle.vertex2, triangle.vertex3, triangle.vertex1));
    m_faces[triangle.vertex3].push_back(Triangle(triangle.vertex3, triangle.vertex1, triangle.vertex2));
}

Vector3d Surface::getPos(int index)
{
    return m_points[index];
}

void Surface::changePos(int index, Vector3d& newPoint)
{
    m_points[index] = newPoint;
    return;
}

void Surface::curvatures(int index, double& gaussCurvature, double& meanCurvature)
{
    // Pre-initialisations to avoid doing it every loop
    Vector3d edge;
    double length;
    Vector3d triangleEdge1;
    Vector3d triangleEdge2;
    double edge1Length;
    double edge2Length;
    double vertexAngle;
    Vector3d triangleEdge3;
    Vector3d triangleEdge4;
    Vector3d norm1;
    Vector3d norm2;
    // End

    auto key_selector = [](auto pair){return pair.first;};
            double curvatureSum = 0;
            std::vector<int> completedEdges;
            Vector3d currentPoint = getPos(index);
            std::unordered_map<int, std::vector<Triangle>> pairTriangles;
            for (Triangle triangle : getNeighbourTriangles(index)) {
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
                    edge = getPos(trianglePairIndex) - currentPoint;
                    length = edge.norm();
                    triangleEdge1 = getPos(triangle1.vertex2) - currentPoint;
                    triangleEdge2 = getPos(triangle1.vertex3) - currentPoint;
                    edge1Length = (triangleEdge1).norm();
                    edge2Length = (triangleEdge2).norm();
                    vertexAngle = triangleEdge1.normalized().dot(triangleEdge2.normalized());
                    triangleEdge3 = getPos(triangle2.vertex2) - currentPoint;
                    triangleEdge4 = getPos(triangle2.vertex3) - currentPoint;
                    norm1 = crossProd(triangleEdge1, triangleEdge2).normalized();
                    norm2 = crossProd(triangleEdge3, triangleEdge4).normalized();
                    area += 0.125 * edge1Length * edge2Length * std::sqrt(1-std::pow(vertexAngle,2));
                    angleSum += std::acos(vertexAngle); 
                    curvatureSum += length * std::atan2(edge.normalized().dot(crossProd(norm2,norm1)), norm2.dot(norm1));
                }
            }
    meanCurvature = (0.25 * curvatureSum)/area;
    gaussCurvature = (2 * M_PI - angleSum)/area;
    return;
}

std::vector<Triangle> Surface::getNeighbourTriangles(int index)
{
    return m_faces[index];
}

Vector3d Surface::crossProd(Vector3d &a, Vector3d &b)
{
    return Vector3d(a[1]*b[2] - a[2]*b[1] ,a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}
