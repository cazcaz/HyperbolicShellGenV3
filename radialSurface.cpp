#include "radialSurface.h"

RadialSurface::RadialSurface() : m_curveCount(0) , m_iterCount(0)
{
}

RadialSurface::~RadialSurface()
{
}

void RadialSurface::addCurve(std::vector<Vector3d> newCurve)
{

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
    }
}

Vector3d RadialSurface::getPoint(int curve, double s)
{
    double pointLocation =  s*double(m_curveLengths[curve])/(2 * M_PI);
    //s is the radial parameter in [0, 2 * pi)
    int prevIndex = int(pointLocation);
    double lineRatio = pointLocation - prevIndex;
    return getPoint(curve, correctIndex(curve, prevIndex)) + lineRatio * (getPoint(curve, correctIndex(curve, prevIndex + 1)) - getPoint(curve, correctIndex(curve, prevIndex)));
}

Vector3d RadialSurface::getPoint(int curve, int index)
{
    int newIndex = curveStartIndex(curve) + index;
    return getPos(newIndex);
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