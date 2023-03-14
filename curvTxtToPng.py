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
        #Each data will be two lists of the gaussian curvs and mean curvs
        datas = []
        surfaces = []
        for filename in os.listdir(fullInputPath):
            if len(filename) >= 13:
                if filename[-13:] == 'curvature.txt':
                    fileNames.append(filename[:-4])
                    with open(os.path.join(fullInputPath, filename), 'r') as f:
                        data = [[],[]]
                        curves=[]
                        currentFile = f.readlines()
                        curveCurvaturesStrings = currentFile[0].split("|")
                        for curveCurvatureString in curveCurvaturesStrings:
                            curvePairStrings = curveCurvatureString.split(":")
                            curve = []
                            for curvePairString in curvePairStrings:
                                currentData = [[],[]]
                                currentValues = curvePairString.split(" ")
                                currentData[0].append(float(currentValues[0]))
                                currentData[1].append(float(currentValues[1]))
                                data[0].append(float(currentValues[0]))
                                data[1].append(float(currentValues[1]))
                                curve.append(currentData)
                            curves.append(curve)
                        datas.append(data)
                        surfaces.append(curves)

        os.chdir(fullInputPath)

        for surface, fileName in zip(surfaces, fileNames):
            thetaInputs = []
            radiusInputs = []
            currentRadius = 1
            gaussColours = []
            meanColours = []
            allGaussColours = []
            allMeanColours = []
            for curve in surface:
                currentCurveLength = len(curve)
                currentGaussCurvatures = []
                currentMeanCuravtures = []
                thetaChange = 2 * pi/currentCurveLength
                currentTheta = 0
                for data in curve:
                    thetaInputs.append(currentTheta)
                    radiusInputs.append(currentRadius)
                    currentTheta += thetaChange
                    currentGaussCurvatures.append(data[0][0])
                    currentMeanCuravtures.append(data[1][0])
                    allGaussColours.append(data[0][0])
                    allMeanColours.append(data[1][0])
                gaussColours.append(currentGaussCurvatures)
                meanColours.append(currentMeanCuravtures)
                currentRadius += 0.1

            q1g, q3g = np.percentile(allGaussColours, [25, 75])
            q1m, q3m = np.percentile(allMeanColours, [25, 75])
            iqrg = q3g - q1g
            iqrm = q3m - q1m
            threshold = 1.5
            upper_boundg = q3g + threshold*iqrg
            upper_boundm = q3m + threshold*iqrm
            lower_boundg = q1g - threshold*iqrg
            lower_boundm = q1m - threshold*iqrm

            longestLength = len(gaussColours[-1])
            interpolatedGaussColours = []
            interpolatedMeanColours = []
            for gaussColourList in gaussColours:
                thetas = np.linspace(0, 2*np.pi, len(gaussColourList))
                newThetas = np.linspace(0, 2*np.pi, longestLength)
                interpolationFunction = interp1d(thetas, np.array(gaussColourList), kind = 'cubic')
                interpolatedGauss = interpolationFunction(newThetas)
                interpolatedGaussColours.append(interpolatedGauss)
            
            for meanColourList in meanColours:
                thetas = np.linspace(0, 2*np.pi, len(meanColourList))
                newThetas = np.linspace(0, 2*np.pi, longestLength)
                interpolationFunction = interp1d(thetas, np.array(meanColourList), kind = 'cubic')
                interpolatedMean = interpolationFunction(newThetas)
                interpolatedMeanColours.append(interpolatedMean)

            gaussColours = np.clip(np.array(interpolatedGaussColours), lower_boundg, upper_boundg).T
            meanColours = np.clip(np.array(interpolatedMeanColours), lower_boundm, upper_boundm).T

            theta, r = np.mgrid[0:2*np.pi:len(surface[-1])*1j, min(radiusInputs):max(radiusInputs):len(surface) *1j]
            

            fig = plt.figure(figsize=(10,6))
            ax1 = fig.add_subplot(121, projection='polar')
            ax2 = fig.add_subplot(122, projection='polar')
            ax1.grid(False)
            ax1.set_rgrids([])
            ax1.set_thetagrids([])
            ax1.set_rorigin(0.5)
            ax2.grid(False)
            ax2.set_rgrids([])
            ax2.set_thetagrids([])
            ax2.set_rorigin(0.5)
            cmap = plt.cm.get_cmap('bwr')
            sc1 = ax1.pcolormesh(theta, r, gaussColours, cmap=cmap)
            sc2 = ax2.pcolormesh(theta, r, meanColours, cmap=cmap)
            cbar1 = plt.colorbar(sc1, ax=ax1, orientation='horizontal')
            cbar2 = plt.colorbar(sc2, ax=ax2, orientation='horizontal')
            cbar1.set_label('Discrete Gaussian Curvatures')
            cbar2.set_label('Discrete Mean Curvatures')
            ax1.set_title('Gaussian Curvature')
            ax2.set_title('Mean Curvature')
            fig.savefig(os.path.join(fullInputPath, fileName + 'sGraph.png'))
            
        for data, fileName in zip(datas, fileNames):
            dataset1 = np.array(data[0])
            mean1 = np.mean(dataset1)
            q1s1 = np.percentile(dataset1,25)
            q2s1 = np.percentile(dataset1,75)

            dataset2 = np.array(data[1])
            mean2 = np.mean(dataset2)
            q1s2 = np.percentile(dataset2,25)
            q2s2 = np.percentile(dataset2,75)

            fig, axs = plt.subplots(2,1, figsize=(6,8))
            axs[0].boxplot(dataset1, meanprops=dict(color='red'))
            axs[0].set_title('Gaussian Curvature')
            axs[0].set_ylabel('Curvature')
            axs[0].text(0.95,0.95, 'Mean:' + str(round(mean1,4)), transform=axs[0].transAxes, ha='right', va='top')
            axs[0].text(0.95,0.90, 'IQRL:' + str(round(q1s1,4)), transform=axs[0].transAxes, ha='right', va='top')
            axs[0].text(0.95,0.85, 'IQRU:' + str(round(q2s1,4)), transform=axs[0].transAxes, ha='right', va='top')

            axs[1].boxplot(dataset2, meanprops=dict(color='red'))
            axs[1].set_title('Mean Curvature')
            axs[1].set_ylabel('Curvature')
            axs[1].text(0.95,0.95, 'Mean:' + str(round(mean2,4)), transform=axs[1].transAxes, ha='right', va='top')
            axs[1].text(0.95,0.90, 'IQRL:' + str(round(q1s2,4)), transform=axs[1].transAxes, ha='right', va='top')
            axs[1].text(0.95,0.85, 'IQRU:' + str(round(q2s2,4)), transform=axs[1].transAxes, ha='right', va='top')

            fig.savefig(os.path.join(fullInputPath, fileName + '.png'))
            os.remove(fileName + '.txt')
            count+= 1
        os.chdir(current)

print("Done, " , count , " .txt files converted to .png.")
