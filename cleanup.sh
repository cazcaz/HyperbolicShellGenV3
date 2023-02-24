#!/bin/bash

# Remove all files and directories in OutputSurfaceTxts except for .gitkeep
rm -rf OutputSurfaceTxts/* ! *.gitkeep

# Remove all files and directories in OutputSurfaceMeshes except for .gitkeep
rm -rf OutputSurfaceMeshes/* ! *.gitkeep

# Remove all files and directories in OutputCurvatures except for .gitkeep
rm -rf CurvatureTxts/* ! *.gitkeep

# Remove all files and directories in OutputCurvaturesPngs except for .gitkeep
rm -rf CurvatureOutputPngs/* ! *.gitkeep

