#pragma once

struct ShellParams {
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 1;
    int resolution = 100;
    int expansions = 50;
    double extensionLength = 0.1;
    double meanStiffness = 0.01;
    double gaussStiffness = 0.1;
    double desiredCurvature = -0.01;
    double strainCoeff = 10000;
    double lengthStiffness = 1;
    double bendingStiffness = 0.001;
    double period = 0;
    bool saveEveryFrame = false;
};
