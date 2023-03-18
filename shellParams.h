#pragma once
#include <iostream>

struct ShellParams
{
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 0.33;
    int resolution = 50;
    int expansions = 100;
    double extensionLength = 0.05;
    double desiredCurvature = -0.01;
    double springCoeff = 10000000;
    double lengthStiffness = 0.1;
    double bendingStiffness = 10000;
    double sharpBend = 0;
    double period = 5;
    bool saveEveryFrame = false;
    int surfaceIndex = 0;
};
