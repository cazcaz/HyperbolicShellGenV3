#pragma once

struct ShellParams {
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius = 1;
    int resolution = 20;
    int expansions = 50;
    double extensionLength = 0.1;
    double meanStiffness = 0.0001;
    double gaussStiffness = 0.1;
    double desiredCurvature = -0.01;
    double strainCoeff = 100000;
    double lengthStiffness = 10;
};
