#pragma once

struct ShellParams {
  double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    int resolution = 50;
    int expansions = 10;
    double extensionLength = 0.1;
    double stiffnessRatio= 0.1;
    double desiredCurvature = -0.000001;
};
