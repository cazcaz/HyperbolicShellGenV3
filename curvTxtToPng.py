import os
import numpy as np
import matplotlib.pyplot as plt

current = os.getcwd()
searchPath = os.path.join(current, "CurvatureTxts")
count = 0
for folder in os.listdir(searchPath):
    fullInputPath = os.path.join(searchPath, folder)
    if os.path.isdir(fullInputPath):
        fullOutputPath = os.path.join(current, "CurvatureOutputPngs", folder)
        fileNames = []
        #Each data will be two lists of the gaussian curvs and mean curvs
        datas = []
        for filename in os.listdir(fullInputPath):
            if not filename.endswith('.txt'):
                continue
            fileNames.append(filename[:-4])
            with open(os.path.join(fullInputPath, filename), 'r') as f:
                data = [[],[]]
                currentFile = f.readlines()
                curvaturePairs = currentFile[0].split(":")
                for pair in curvaturePairs:
                    currentValues = pair.split(" ")
                    data[0].append(float(currentValues[0]))
                    data[1].append(float(currentValues[1]))
                datas.append(data)
        os.chdir(fullOutputPath)
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

            fig.savefig(os.path.join(fullOutputPath, fileName + '.png'))
            count+= 1
        os.chdir(current)

print("Done, " , count , " .txt files converted to .png.")
