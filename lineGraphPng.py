import os
import numpy as np
import matplotlib.pyplot as plt

current = os.getcwd()
searchPath = os.path.join(current, "Surfaces")
count = 0
persistData = []
for folder in os.listdir(searchPath):
    fullInputPath = os.path.join(searchPath, folder)
    if os.path.isdir(fullInputPath):
        fileNames = []
        #Each data will be two lists of the gaussian curvs and mean curvs
        datas = []
        parameters = []
        for filename in os.listdir(fullInputPath):
            if len(filename) >= 12:
                if filename[-12:] == 'linePlot.txt':
                    fileNames.append(filename[:-12])
                    with open(os.path.join(fullInputPath, filename), 'r') as f:
                        currentFile = f.readlines()
                        paramterStringList = []
                        for line in currentFile:
                            if line[0] != "?":
                                paramterStringList.append(line[:-1])
                        parameters.append(paramterStringList)
                        curvaturePairs = currentFile[-1][1:].split(":")
                        data = [[] for _ in range(len(curvaturePairs[0].split(" ")[0]))]
                        for pair in curvaturePairs:
                            currentValues = pair.split(" ")
                            for dataIndex in range(len(currentValues)):
                                data[dataIndex].append(float(currentValues[dataIndex]))
                        datas.append(data)
        os.chdir(fullInputPath)
        for parameterStrings, data, fileName in zip(parameters, datas, fileNames):
            xvalues = np.array(data[0])
            yvalues = [np.array(ydata) for ydata in data[1:]]
            if fileName == "energyProfile":
                persistData.append([parameterStrings[3], xvalues, yvalues[0]])
            fig, ax = plt.subplots()
            ax.set_title(parameterStrings[0])
            axisLabels = parameterStrings[1].split("~")
            ax.set_ylabel(axisLabels[1])
            ax.set_xlabel(axisLabels[0])
            currentYTextCoord = 1.15
            for parameter in parameterStrings[4:]:
                fig.text(-0.15,currentYTextCoord, parameter, transform=ax.transAxes, ha='left', va='top',  fontsize=8)
                currentYTextCoord -= 0.025
            legends = parameterStrings[2].split("~")
            for yvalueset, legend in zip(yvalues, legends):
                ax.plot(xvalues,yvalueset, label = legend)

            ax.legend(loc='upper left')
            fig.savefig(os.path.join(fullInputPath, fileName + '.png'))
            os.remove(fileName + 'linePlot.txt')
            count+= 1
        os.chdir(current)
if len(persistData) != 0:
    currentPlottingCount = 1
    totalPlottingCount = 0
    while totalPlottingCount < len(persistData):
        fig, ax = plt.subplots()
        ax.set_title("Combined Energy Plots")
        ax.set_ylabel("Energy")
        ax.set_xlabel("Radius")
        endCount = totalPlottingCount + 5
        if totalPlottingCount+5 > len(persistData)-1:
            endCount = len(persistData)-1
        for dataTriples in persistData[totalPlottingCount:endCount]:
            ax.plot(dataTriples[1], dataTriples[2], label = dataTriples[0])
        ax.legend(loc='upper left')
        fig.savefig(os.path.join(searchPath, 'combinedEnergy' +str(currentPlottingCount)+ '.png'))
        currentPlottingCount += 1
        totalPlottingCount += 5
    count += 1
print("Done, " , count , " .txt files converted to .png.")
