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
    m_faces[triangle.vertex2].push_back(triangle);
    m_faces[triangle.vertex3].push_back(triangle);


}

Vector3d Surface::getPos(int index)
{
    return m_points[index];
}

double Surface::meanCurvature(int index)
{
    std::vector<double> edgeLengths;
    std::vector<int> completedEdges;
    Vector3d currentPoint = getPos(index);
    int uniqueNeighbour1;
    int uniqueNeighbour2;
    for (Triangle triangle : m_faces[index]) {
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
        bool unique1 = true;
        bool unique2 = true;
        for (int completedEdge : completedEdges) {
            if (completedEdge == uniqueNeighbour1) {
                unique1 = false;
            } 
            if (completedEdge == uniqueNeighbour2) {
                unique2 = false;
            }
        }
        if (unique1) {
            edgeLengths.push_back((currentPoint - getPos(uniqueNeighbour1)).norm());
            completedEdges.push_back(uniqueNeighbour1);
            std::cout << uniqueNeighbour1 << std::endl;
        }
        if (unique2) {
            edgeLengths.push_back((currentPoint - getPos(uniqueNeighbour2)).norm());
            completedEdges.push_back(uniqueNeighbour2);
            std::cout << uniqueNeighbour2 << std::endl;
        }

    }
    std::cout << edgeLengths.size() << std::endl;
    return 0;
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
        
        edge1 = (vertex1 - midPoint).normalized();
        edge2 = (vertex2 - midPoint).normalized();
        angleSum += std::acos(edge1.dot(edge2));
    }
    
    return 2 * M_PI - angleSum;
}
