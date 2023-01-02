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

double Surface::meanCurvature(int index)
{
    auto key_selector = [](auto pair){return pair.first;};
    double curvatureSum = 0;
    std::vector<int> completedEdges;
    Vector3d currentPoint = getPos(index);
    std::unordered_map<int, std::vector<Triangle>> pairTriangles; 
    int uniqueNeighbour1;
    int uniqueNeighbour2;
    for (Triangle triangle : m_faces[index]) {
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
    for (int trianglePairIndex : keys) {
        Triangle triangle1 = pairTriangles[trianglePairIndex][0];
        Triangle triangle2 = pairTriangles[trianglePairIndex][1];
        if (triangle1.vertex2 == triangle2.vertex3) {
            triangle1 = triangle2;
            triangle2 = pairTriangles[trianglePairIndex][0];
        }
        Vector3d edge = getPos(trianglePairIndex) - currentPoint;
        double length = edge.norm();
        Vector3d triangleEdge1 = getPos(triangle1.vertex2) - currentPoint;
        Vector3d triangleEdge2 = getPos(triangle1.vertex3) - currentPoint;
        double edge1Length = (triangleEdge1).norm();
        double edge2Length = (triangleEdge2).norm();
        double vertexAngle = std::acos(triangleEdge1.normalized().dot(triangleEdge2.normalized()));
        Vector3d triangleEdge3 = getPos(triangle2.vertex2) - currentPoint;
        Vector3d triangleEdge4 = getPos(triangle2.vertex3) - currentPoint;
        Vector3d norm1 = crossProd(triangleEdge1, triangleEdge2).normalized();
        Vector3d norm2 = crossProd(triangleEdge3, triangleEdge4).normalized();
        area += 0.125 * edge1Length * edge2Length * std::sin(vertexAngle);
        curvatureSum += length * std::atan2(edge.normalized().dot(crossProd(norm2,norm1)), norm2.dot(norm1));//1/triangleArea * length * std::atan2(edge.normalized().dot(crossProd(norm2,norm1)), norm2.dot(norm1));
    }
    return (0.25 * curvatureSum)/area;
}

double Surface::gaussCurvature(int index)
{
    double angleSum = 0;
    std::vector<Triangle> pointTriangles = m_faces[index];

    int uniqueNeighbour1;
    int uniqueNeighbour2;
    Vector3d vertex1;
    Vector3d vertex2;
    Vector3d midPoint = getPos(index);
    Vector3d edge1;
    Vector3d edge2;
    double area = 0;
    for (Triangle triangle : pointTriangles) {
        if (triangle.vertex1 == index) {
            uniqueNeighbour1 = triangle.vertex2;
            uniqueNeighbour2 = triangle.vertex3;
        } else if (triangle.vertex2 == index) {
            uniqueNeighbour1 = triangle.vertex1;
            uniqueNeighbour2 = triangle.vertex3;
        } else {
            uniqueNeighbour1 = triangle.vertex1;
            uniqueNeighbour2 = triangle.vertex2;
        }
        vertex1 = getPos(uniqueNeighbour1);
        vertex2 = getPos(uniqueNeighbour2);
        
        double edge1Length = (vertex1 - midPoint).norm();
        double edge2Length = (vertex2 - midPoint).norm();
        edge1 = (vertex1 - midPoint).normalized();
        edge2 = (vertex2 - midPoint).normalized();
        double vertexAngle = std::acos(edge1.dot(edge2));
        area += 0.125 * edge1Length * edge2Length * std::sin(vertexAngle);
        angleSum +=  vertexAngle;
    }
    return (2 * M_PI - angleSum)/area;
}


Vector3d Surface::crossProd(Vector3d &a, Vector3d &b)
{
    return Vector3d(a[1]*b[2] - a[2]*b[1] ,a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}
