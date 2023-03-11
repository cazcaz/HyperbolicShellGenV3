#pragma once

struct ShellParams {
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 1;
    int resolution = 50;
    int expansions = 30;
    double extensionLength = 0.1;
    double meanStiffness = 0.01;
    double gaussStiffness = 0.1;
    double desiredCurvature = -0.01;
    double strainCoeff = 1000000;
    double lengthStiffness = 1;
    double bendingStiffness = 10;
    double period = 0;
    bool saveEveryFrame = false;
};
