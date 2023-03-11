#pragma once

#include "surface.h"

class RadialSurface : public Surface {
    public:
        RadialSurface();
        ~RadialSurface();
        void addCurve(std::vector<Vector3d> newCurve);
        Vector3d getPoint(int curve, double s);
        Vector3d getPoint(int curve, int index);
        Vector3d getEdgeNormal(int curve, double s);
        Vector3d getEdgeTangent(int curve, double s);
        std::vector<Vector3d>getCurve(int curveNumber);
        void changeCurvePoint(int curve, int index, Vector3d& newPoint);
        int getCurveCount();
        int getIterCount();
        void increaseIterCount();
        void addIterCount(int count);
        int curveStartIndex(int curveNum);
        int correctIndex(int curve, int index);
        int getCurveSize(int curveNum);
        double getCurveLength(int curveNum);
        
    private:
        int m_curveCount;
        std::vector<int> m_curveLengths;
        std::vector<std::vector<double>> m_segmentLengths;
        int m_iterCount;
};