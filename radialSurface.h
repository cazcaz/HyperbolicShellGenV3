#pragma once

#include "surface.h"

class RadialSurface : public Surface {
    public:
        RadialSurface();
        ~RadialSurface();
        void addCurve(std::vector<Vector3d> newCurve);
        Vector3d getPoint(int curve, double s);
        Vector3d getPoint(int curve, int index);
        std::vector<Vector3d>getCurve(int curveNumber);
        int getCurveCount();
    private:
        int curveStartIndex(int curveNum);
        int correctIndex(int curve, int index);
        int m_curveCount;
        std::vector<int> m_curveLengths;
};