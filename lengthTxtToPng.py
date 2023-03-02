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
        for filename in os.listdir(fullInputPath):
            if len(filename) >= 17:
                if filename[-17:] == 'lengthProfile.txt':
                    fileNames.append(filename[:-4])
                    with open(os.path.join(fullInputPath, filename), 'r') as f:
                        data = [[],[],[]]
                        currentFile = f.readlines()
                        curvaturePairs = currentFile[0].split(":")
                        for pair in curvaturePairs:
                            currentValues = pair.split(" ")
                            data[0].append(float(currentValues[0]))
                            data[1].append(float(currentValues[1]))
                            data[2].append(float(currentValues[2]))
                        datas.append(data)
        os.chdir(fullInputPath)
        for data, fileName in zip(datas, fileNames):
            xvalues = np.array(data[0])
            curveRealLengths = np.array(data[1])
            curveDesiredLengths = np.array(data[2])

            fig, ax = plt.subplots()
            ax.plot(xvalues,curveRealLengths, label = "Actual Radius Lengths")
            ax.plot(xvalues,curveDesiredLengths, '--', label = "Desired Radius Lengths")
            ax.set_title('Curve Lengths')
            ax.set_ylabel('Length')
            ax.set_xlabel('Radius')
            ax.legend(loc='upper left')
            fig.savefig(os.path.join(fullInputPath, fileName + '.png'))
            os.remove(fileName + '.txt')
            count+= 1
        os.chdir(current)

print("Done, " , count , " .txt files converted to .png.")
