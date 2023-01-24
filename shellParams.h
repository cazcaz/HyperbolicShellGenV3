#pragma once

struct ShellParams {
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 1;
    int resolution = 100;
    int expansions = 50;
    double extensionLength = 0.1;
    double stiffnessRatio= 10;
    double desiredCurvature = 0.1;
    double strainCoeff = 1000;
    double lengthStiffness = 100;
};
