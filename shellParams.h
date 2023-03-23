#pragma once
#include <iostream>

struct ShellParams
{
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 0.33;
    double bendingStiffness = 5000;
    double desiredCurvature = -2.5;
    int resolution = 50;
    int expansions = 60;
    double extensionLength = 0.05;
    double springCoeff = 10000000;
    double sharpBend = 0;
    double period = 5;
    double lengthStiffness = 0.1;
    bool saveEveryFrame = false;
    int surfaceIndex = 0;
};
