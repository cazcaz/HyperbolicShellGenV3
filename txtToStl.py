from stl import mesh
import os
import numpy as np

current = os.getcwd()
searchPath = os.path.join(current, "Surfaces")
count = 0
for folder in os.listdir(searchPath):
    fullInputPath = os.path.join(searchPath, folder)
    if folder == ".gitkeep":
        continue
    os.chdir(fullInputPath)
    if os.path.isdir(fullInputPath):
        surfaces = []
        fileNames = []
        for filename in os.listdir(fullInputPath):
            if len(filename) >= 11:
                if filename[-11:] == 'surface.txt':
                    fileNames.append(filename[:-4])
                    with open(os.path.join(fullInputPath, filename), 'r') as f:
                        currentSurface = []
                        currentFile = f.readlines()
                        curvesStrs = currentFile[0].split("|")
                        for curveString in curvesStrs:
                            pointData = curveString.split(":")
                            pointString = pointData[0].split(",")
                            triangleStrings = pointData[1].split("~")
                            point = [float(coord) for coord in pointString]
                            neighbours = [[int(index) for index in indices.split(",")] for indices in triangleStrings]
                            currentSurface.append([point, neighbours])
                        listVertices = []
                        listTriangles = []
                        point = 0
                        for pointTrianglePair in currentSurface:
                            listVertices.append(pointTrianglePair[0])
                            for triangle in pointTrianglePair[1]:
                                validTriangle = True
                                for index in triangle:
                                    if point > index:
                                        validTriangle = False
                                        break
                                if validTriangle:
                                    listTriangles.append(triangle)
                            point += 1
                        vertices = np.array(listVertices)
                        triangles = np.array(listTriangles)
                        surfaceMesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
                        for i, f in enumerate(triangles):
                            for j in range(3):
                                surfaceMesh.vectors[i][j] = vertices[f[j],:]
                        surfaceMesh.save(os.path.abspath(filename[:-4] + '.stl'))
                        os.remove(filename)
                        count+= 1

print("Done, " , count , " .txt files converted to .stl.")