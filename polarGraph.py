import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
        radiusDatas = []
        parameters = []
        for filename in os.listdir(fullInputPath):
            if len(filename) >= 13:
                if filename[-13:] == 'polarPlot.txt':
                    fileNames.append(filename[:-13])
                    with open(os.path.join(fullInputPath, filename), 'r') as f:
                        currentFile = f.readlines()
                        paramterStringList = []
                        for line in currentFile:
                            if line[0] != "?":
                                paramterStringList.append(line[:-1])
                        parameters.append(paramterStringList)
                        curves=[]
                        radiusData=[]
                        curveDataStrings = currentFile[-1][1:].split("|")
                        for curveDataString in curveDataStrings:
                            IndDispString = curveDataString.split(" ")[1].split(",")
                            radiusData.append(float(curveDataString.split(" ")[0]))
                            curve = []
                            for IndDisp in IndDispString:
                                curve.append(float(IndDisp))
                            curves.append(curve)
                        radiusDatas.append(radiusData)
                        surfaces.append(curves)

        os.chdir(fullInputPath)

        for parameterStrings, surface, radiusData, fileName in zip(parameters, surfaces, radiusDatas, fileNames):
            thetaInputs = []
            polarData = []
            allPolarData = []
            for curve in surface:
                currentPolarData = []
                currentCurveLength = len(curve)
                thetaChange = 2 * pi/currentCurveLength
                currentTheta = 0
                for data in curve:
                    thetaInputs.append(currentTheta)
                    currentTheta += thetaChange
                    currentPolarData.append(data)
                    allPolarData.append(data)
                polarData.append(currentPolarData)

            longestLength = len(polarData[-1])
            interpolatedPolarDatas = []
            for polarDataList in polarData:
                thetas = np.linspace(0, 2*np.pi, len(polarDataList))
                newThetas = np.linspace(0, 2*np.pi, longestLength)
                interpolationFunction = interp1d(thetas, np.array(polarDataList), kind = 'cubic')
                interpolatedPolarData = interpolationFunction(newThetas)
                print(np.mean(interpolatedPolarData))
                interpolatedPolarDatas.append(interpolatedPolarData)

            if parameterStrings[3] == "1":
                q1, q3 = np.percentile(allPolarData, [25,75])
                iqr = q3 - q1
                threshold = 1.5
                upper_bound = q3 + threshold*iqr
                lower_bound = q1 - threshold*iqr
                interpolatedPolarDatas = np.clip(np.array(interpolatedPolarDatas), lower_bound, upper_bound)
            
            thetaValuesFinal = np.linspace(0, 2*np.pi, len(surface[-1]), endpoint = False)
            theta, r = np.meshgrid(thetaValuesFinal, radiusData)
            z = np.array(interpolatedPolarDatas)
            print(parameterStrings, np.mean(z))
            fig = plt.figure(figsize=(5,6))
            ax1 = fig.add_subplot(111, projection='polar')
            ax1.grid(False)
            ax1.set_rgrids([])
            ax1.set_thetagrids([])
            ax1.set_rorigin(radiusData[0])
            currentYTextCoord = 1.2
            for parameter in parameterStrings[4:]:
                fig.text(-0.25,currentYTextCoord, parameter, transform=ax1.transAxes, ha='left', va='top')
                currentYTextCoord -= 0.05
            cmap = plt.cm.get_cmap(parameterStrings[2])
            sc1 = ax1.pcolormesh(theta, r, z, norm=colors.CenteredNorm(),cmap=cmap)
            cbar1 = plt.colorbar(sc1, ax=ax1, orientation='horizontal')
            cbar1.set_label(parameterStrings[1])
            ax1.set_title(parameterStrings[0])
            fig.savefig(os.path.join(fullInputPath, fileName + '.png'))
            os.remove(fileName + 'polarPlot.txt')
            count+= 1
        os.chdir(current)

print("Done, " , count , " .txt files converted to .png.")