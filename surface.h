#pragma once

#include <Eigen/core>
#include <vector>
#include <iostream>
#include <unordered_map>

using Eigen::Vector3d;

class Surface {
    public:
        Surface();
        ~Surface();

        void addPoint(Vector3d& newPoint);
        void addNeighbour(int position, int newNeighboursIndex);
        Vector3d getPos(int index);
        friend std::ostream& operator<< (std::ostream& os, Surface& surface){
        for (int i=0; i < surface.m_points.size(); i++){
            Vector3d point = surface.getPos(i);
            os << point[0] << "," << point[1] << "," << point[2] << ":";
            std::vector<int> neighbours = surface.m_neighbours[i];
            for (int j = 0; j < neighbours.size();j++) {
                if (j == neighbours.size()-1) {
                    if (i == surface.m_points.size() - 1){
                        os << neighbours[j];
                    } else {
                        os << neighbours[j] << "|";
                    }
                } else {
                    os << neighbours[j] << ",";
                }
            }
        }
        return os;
        };
    private:
        std::vector<Vector3d> m_points;
        std::unordered_map<int, std::vector<int>> m_neighbours;
        int m_pointCount;
};
