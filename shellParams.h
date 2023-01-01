#pragma once

struct ShellParams {
  double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    int resolution = 50;
    int expansions = 10;
    double extensionLength = 0.05;
    double stiffLengthRatioCircum = 0.1;
    double stiffLengthRatioRadial = 0.1;
    double desiredCurvature = 0.001;
};
