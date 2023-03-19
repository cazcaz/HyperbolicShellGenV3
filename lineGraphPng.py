import os
import numpy as np
import matplotlib.pyplot as plt

current = os.getcwd()
searchPath = os.path.join(current, "Surfaces")
count = 0
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

            fig, ax = plt.subplots()
            ax.set_title(parameterStrings[0])
            axisLabels = parameterStrings[1].split("~")
            ax.set_ylabel(axisLabels[1])
            ax.set_xlabel(axisLabels[0])
            currentYTextCoord = 1.15
            for parameter in parameterStrings[3:]:
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

print("Done, " , count , " .txt files converted to .png.")
