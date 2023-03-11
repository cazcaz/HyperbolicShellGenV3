#include "radialSurface.h"

RadialSurface::RadialSurface() : m_curveCount(0) , m_iterCount(0)
{
}

RadialSurface::~RadialSurface()
{
}

void RadialSurface::addCurve(std::vector<Vector3d> newCurve)
{
    std::vector<double> segmentLengths;
    for (Vector3d point : newCurve) {
        addPoint(point);
    }
    bool initCurve = false;
    if(m_curveLengths.size() == 0){
        initCurve = true;
    }
    int newCurveSize = newCurve.size();
    m_curveLengths.push_back(newCurveSize);
    m_curveCount += 1;
    if (!initCurve) {
        int currentConnectionIndex = curveStartIndex(m_curveCount-2);
        double pointGapOld = 1/double(m_curveLengths[m_curveCount-2]);
        double pointGapNew = 1/double(newCurveSize);
        double rangeCounterOld = 0.5* pointGapOld;
        double rangeCounterNew = 0;
        for(int i = 0; i < newCurveSize; i++) {
            if (i < newCurveSize-1){
                segmentLengths.push_back((newCurve[i] - newCurve[i+1]).norm());
            } else {
                segmentLengths.push_back((newCurve[i] - newCurve[0]).norm());
            }
            int currentIndex = curveStartIndex(m_curveCount-1) + i;
            if (rangeCounterNew > rangeCounterOld){
                rangeCounterOld += pointGapOld;
                currentConnectionIndex += 1;
                addTriangle(Triangle(currentIndex, currentConnectionIndex - 1, currentConnectionIndex));
            }
            addTriangle(Triangle(currentIndex, currentConnectionIndex ,curveStartIndex(m_curveCount - 1) +  correctIndex(m_curveCount-1, currentIndex+1 - curveStartIndex(m_curveCount - 1))));
            rangeCounterNew += pointGapNew;
        }
        addTriangle(Triangle(correctIndex(m_curveCount - 2, -1) + curveStartIndex(m_curveCount-2),  curveStartIndex(m_curveCount-2),curveStartIndex(m_curveCount - 1)));
    } else {
        for(int i = 0; i < newCurveSize; i++) {
            if (i < newCurveSize-1){
                segmentLengths.push_back((newCurve[i] - newCurve[i+1]).norm());
            } else {
                segmentLengths.push_back((newCurve[i] - newCurve[0]).norm());
            }
        }
    }
    m_segmentLengths.push_back(segmentLengths);
}

Vector3d RadialSurface::getPoint(int curve, double s)
{
    double curveLength = getCurveLength(curve);
    int curvePointCount = getCurveSize(curve);
    double rescaledS = s/(2*M_PI) * curveLength;
    int segmentCount = 0;
    double currentLength = m_segmentLengths[curve][0];
    while (currentLength < rescaledS) {
        segmentCount += 1;
        currentLength += m_segmentLengths[curve][segmentCount];
    }
    currentLength -= m_segmentLengths[curve][segmentCount];
    double tParam = rescaledS - currentLength;
    int secondIndex = segmentCount+1;
    if (secondIndex >= curvePointCount) {
        secondIndex -= curvePointCount;
    }
    Vector3d edgeStart = getPoint(curve, segmentCount);
    Vector3d edgeEnd = getPoint(curve, secondIndex);
    Vector3d finalPoint = edgeStart + tParam/m_segmentLengths[curve][segmentCount] * (edgeEnd - edgeStart);
    return finalPoint;
}

Vector3d RadialSurface::getPoint(int curve, int index)
{
    int newIndex = curveStartIndex(curve) + index;
    return getPos(newIndex);
}

Vector3d RadialSurface::getEdgeNormal(int curve, double s)
{
    double curveLength = getCurveLength(curve);
    int curvePointCount = getCurveSize(curve);
    double rescaledS = s/(2*M_PI) * curveLength;
    int segmentCount = 0;
    double currentLength = m_segmentLengths[curve][0];
    while (currentLength < rescaledS) {
        segmentCount += 1;
        currentLength += m_segmentLengths[curve][segmentCount];
    }
    int secondIndex = segmentCount+1;
    if (secondIndex >= curvePointCount) {
        secondIndex -= curvePointCount;
    }
    std::vector<Triangle> firstsegmentTringles = getNeighbourTriangles(curveStartIndex(curve)+segmentCount);
    Vector3d normal;
    for (Triangle triangle : firstsegmentTringles) {
        if (triangle.vertex3 == secondIndex + curveStartIndex(curve)) {
            Vector3d edge1 = getPos(triangle.vertex3) - getPos(triangle.vertex1);
            Vector3d edge2 = getPos(triangle.vertex2) - getPos(triangle.vertex1);
            normal = Vector3d(edge1[1]*edge2[2] - edge1[2]*edge2[1] ,edge1[2]*edge2[0] - edge1[0]*edge2[2], edge1[0]*edge2[1] - edge1[1]*edge2[0]);
        }
    }
    normal.normalize();

    return normal;
}

Vector3d RadialSurface::getEdgeTangent(int curve, double s)
{
    double curveLength = getCurveLength(curve);
    int curvePointCount = getCurveSize(curve);
    double rescaledS = s/(2*M_PI) * curveLength;
    int segmentCount = 0;
    double currentLength = m_segmentLengths[curve][0];
    while (currentLength < rescaledS) {
        segmentCount += 1;
        currentLength += m_segmentLengths[curve][segmentCount];
    }
    int secondIndex = segmentCount+1;
    if (secondIndex >= curvePointCount) {
        secondIndex -= curvePointCount;
    }
    Vector3d edgeStart = getPoint(curve, segmentCount);
    Vector3d edgeEnd = getPoint(curve, secondIndex);
    Vector3d tangent = (edgeEnd - edgeStart).normalized();
    return tangent;
    
}

std::vector<Vector3d> RadialSurface::getCurve(int curveNumber)
{
    std::vector<Vector3d> curve;
    for (int i = 0; i < m_curveLengths[curveNumber]; i++) {
        curve.push_back(getPoint(curveNumber, i));
    }
    return curve;
}

void RadialSurface::changeCurvePoint(int curve, int index, Vector3d &newPoint)
{
    int curveIndex = curveStartIndex(curve);
    changePos(curveIndex + index, newPoint);
    return;
}

int RadialSurface::getCurveCount()
{
    return m_curveCount;
}

int RadialSurface::getIterCount()
{
    return m_iterCount;
}

void RadialSurface::increaseIterCount()
{
    m_iterCount+=1;
}

void RadialSurface::addIterCount(int count)
{
    m_iterCount += count;
}

int RadialSurface::curveStartIndex(int curveNum)
{
    int curveStartIndex = 0;
    for (int i=0; i < curveNum; i ++) {
        curveStartIndex += m_curveLengths[i];
    }
    return curveStartIndex;
}

int RadialSurface::correctIndex(int curve, int index)
{
    if (index >= m_curveLengths[curve]) {
        return correctIndex(curve, index - m_curveLengths[curve]);
    } else if (index < 0) {
        return correctIndex(curve, index + m_curveLengths[curve]);
    }
    return index;
}

int RadialSurface::getCurveSize(int curveNum)
{
    return m_curveLengths[curveNum];
}

double RadialSurface::getCurveLength(int curveNum) {
    std::vector<Vector3d> curve = getCurve(curveNum);
    double length=0;
    for (int index=0; index < curve.size(); index++) {
        Vector3d point1 = curve[index];
        Vector3d point2;
        if (index == curve.size() - 1) {
            point2 = curve[0];
        } else {
            point2 = curve[index + 1];
        }
        length += (point1 - point2).norm();
    }
    return length;
}