#!/bin/bash

# Remove all files and directories in Surfaces except for .gitkeep
rm -rf Surfaces/* ! *.gitkeep

# Remove all files and directories in OutputSurfaceMeshes except for .gitkeep
rm -rf OutputSurfaceMeshes/* ! *.gitkeep

# Remove all files and directories in OutputCurvatures except for .gitkeep
rm -rf CurvatureTxts/* ! *.gitkeep

# Remove all files and directories in OutputCurvaturesPngs except for .gitkeep
rm -rf CurvatureOutputPngs/* ! *.gitkeep

