import os
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore", message="The get_cmap function was deprecated in Matplotlib 3.7 and will be removed two minor releases later")

current = os.getcwd()
searchPath = os.path.join(current, "Surfaces")
count = 0
for folder in os.listdir(searchPath):
    fullInputPath = os.path.join(searchPath, folder)
    if os.path.isdir(fullInputPath):
        fileNames = []
        surfaces = []
        parameters = []
        for filename in os.listdir(fullInputPath):
            if len(filename) >= 17:
                if filename[-17:] == 'displacements.txt':
                    fileNames.append(filename[:-4])
                    with open(os.path.join(fullInputPath, filename), 'r') as f:
                        currentFile = f.readlines()
                        paramterStringList = []
                        for line in currentFile:
                            if line[0] != "?":
                                paramterStringList.append(line[:-1])
                        parameters.append(paramterStringList)
                        curves=[]
                        curveDispStrings = currentFile[-1][1:].split("|")
                        for curveDispString in curveDispStrings:
                            IndDispString = curveDispString.split(",")
                            curve = []
                            for IndDisp in IndDispString:
                                curve.append(float(IndDisp))
                            curves.append(curve)
                        surfaces.append(curves)

        os.chdir(fullInputPath)

        for parameterStrings, surface, fileName in zip(parameters, surfaces, fileNames):
            thetaInputs = []
            radiusInputs = []
            currentRadius = 1
            displacements = []
            for curve in surface:
                currentDisplacements =[]
                currentCurveLength = len(curve)
                thetaChange = 2 * pi/currentCurveLength
                currentTheta = 0
                for displacement in curve:
                    thetaInputs.append(currentTheta)
                    radiusInputs.append(currentRadius)
                    currentTheta += thetaChange
                    currentDisplacements.append(displacement)
                currentRadius += 0.1
                displacements.append(currentDisplacements)

            longestLength = len(displacements[-1])
            interpolatedDisplacements = []
            for displacementList in displacements:
                thetas = np.linspace(0, 2*np.pi, len(displacementList))
                newThetas = np.linspace(0, 2*np.pi, longestLength)
                interpolationFunction = interp1d(thetas, np.array(displacementList), kind = 'cubic')
                interpolatedDisplacement = interpolationFunction(newThetas)
                interpolatedDisplacements.append(interpolatedDisplacement)


            # q1, q3 = np.percentile(gaussColours, [25, 75])
            # iqr = q3 - q1
            # threshold = 1.5
            # upper_boundg = q3 + threshold*iqr
            # lower_boundg = q1 - threshold*iqr

            # displacements = np.clip(np.array(displacements), lower_bound, upper_bound)
            theta, r = np.mgrid[0:2*np.pi:len(surface[-1])*1j, min(radiusInputs):max(radiusInputs):len(surface) *1j]
            z = np.array(interpolatedDisplacements).T
            fig = plt.figure(figsize=(5,6))
            ax1 = fig.add_subplot(111, projection='polar')
            ax1.grid(False)
            ax1.set_rgrids([])
            ax1.set_thetagrids([])
            ax1.set_rorigin(0.5)
            currentYTextCoord = 1.2
            for parameter in parameterStrings:
                fig.text(-0.25,currentYTextCoord, parameter, transform=ax1.transAxes, ha='left', va='top')
                currentYTextCoord -= 0.05
            cmap = plt.cm.get_cmap('PiYG')
            sc1 = ax1.pcolormesh(theta, r, z, cmap=cmap)
            cbar1 = plt.colorbar(sc1, ax=ax1, orientation='horizontal')
            cbar1.set_label('Input Values')
            ax1.set_title('Input Values')
            fig.savefig(os.path.join(fullInputPath, fileName + '.png'))
            os.remove(fileName + '.txt')
            count+= 1
        os.chdir(current)

print("Done, " , count , " .txt files converted to .png.")
