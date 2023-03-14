#pragma once

struct ShellParams {
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 0.33;
    int resolution = 100;
    int expansions = 100;
    double extensionLength = 0.1;
    double desiredCurvature = -0.01;
    double strainCoeff = 10000;
    double lengthStiffness = 1;
    double bendingStiffness = 10;
    double period = 4;
    bool saveEveryFrame = true;
};
