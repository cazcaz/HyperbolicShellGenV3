#pragma once
#include <iostream>

struct ShellParams
{

    // Parameters that directly effect the surface generation
    double centreX = 0;
    double centreY = 0;
    double centreZ = 0;
    double radius =0.2;
    double bendingStiffness = 1000;
    double radialStiffness = 1;
    double springCoeff = 10000000;
    double desiredCurvature = -3;
    int resolution = 50;
    int expansions = 50;
    double extensionLength = 0.05;
    int surfaceIndex = 0;
    
    // Unused energy function parameters
    double sharpBend = 0;
    double lengthStiffness = 0.1;


    // Enable to save the surface at every point of minimsation, used to make animations of the surface growth
    bool saveEveryFrame = false;

    // Ignore, used and attached to surface parameters for persistance reasons
    bool offsetNeeded = false;
};
