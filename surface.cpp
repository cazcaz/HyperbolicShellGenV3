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
    m_neighbours[m_pointCount-1];
    return;
}

void Surface::addNeighbour(int position, int newNeighboursIndex)
{
    m_neighbours[position].push_back(newNeighboursIndex);
    m_neighbours[newNeighboursIndex].push_back(position);
}

Vector3d Surface::getPos(int index)
{
    return m_points[index];
}

