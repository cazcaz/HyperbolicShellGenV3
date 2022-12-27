#pragma once

#include <eigen3/Eigen/Core>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cmath>

using Eigen::Vector3d;

struct Triangle {
    Triangle(int vert1, int vert2, int vert3) : vertex1{vert1} , vertex2{vert2}, vertex3{vert3} {}
    int vertex1;
    int vertex2;
    int vertex3;
};


class Surface {
    public:
        Surface();
        ~Surface();

        void addPoint(Vector3d& newPoint);
        void addTriangle(Triangle triangle);
        Vector3d getPos(int index);
        double meanCurvature(int index);
        double gaussCurvature(int index);
        friend std::ostream& operator<< (std::ostream& os, Surface& surface){
        for (int i=0; i < surface.m_points.size(); i++){
            Vector3d point = surface.getPos(i);
            os << point[0] << "," << point[1] << "," << point[2] << ":";
            std::vector<Triangle> neighbours = surface.m_faces[i];
            for (int j = 0; j < neighbours.size();j++) {
                if (j == neighbours.size()-1) {
                    if (i == surface.m_points.size() - 1){
                        os << neighbours[j].vertex1 << "," << neighbours[j].vertex2 << "," << neighbours[j].vertex3;
                    } else {
                        os << neighbours[j].vertex1 << "," << neighbours[j].vertex2 << "," << neighbours[j].vertex3;
                        os << "|";
                    }
                } else {
                    os << neighbours[j].vertex1 << "," << neighbours[j].vertex2 << "," << neighbours[j].vertex3;
                    os << "~";
                }
            }
        }
        return os;
        };
    private:
        std::vector<Vector3d> m_points;
        std::unordered_map<int, std::vector<Triangle>> m_faces;
        int m_pointCount;
};
