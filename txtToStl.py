from stl import mesh
import os
import numpy as np

surfaces = []
fileNames = []
current = os.getcwd()
count = 0
os.chdir("..")
os.chdir(current + "/OutputSurfaceTxts")
for filename in os.listdir(os.getcwd()):
    if not filename.endswith('.txt'):
        continue
    fileNames.append(filename[:-4])
    with open(os.path.join(os.getcwd(), filename), 'r') as f:
        currentSurface = []
        currentFile = f.readlines()
        curvesStrs = currentFile[0].split("|")
        for curveString in curvesStrs:
            pointData = curveString.split(":")
            pointString = pointData[0].split(",")
            neighboursString = pointData[1].split(",")
            point = [float(coord) for coord in pointString]
            neighbours = [int(index) for index in neighboursString]
            currentSurface.append([point, neighbours])
        surfaces.append(currentSurface)
os.chdir("..")
current = os.getcwd()
os.chdir(current + "/OutputSurfaceMeshes")
for surface, fileName in zip(surfaces, fileNames):
    listVertices = []
    listTriangles = []
    for i in range(len(surface)):
        listVertices.append(surface[i][0])
        for neighbour in surface[i][1]:
            if neighbour > i:
                for secondNeighbour in surface[i][1]:
                    if secondNeighbour > i and secondNeighbour > neighbour:
                        if secondNeighbour in surface[neighbour][1]:
                            listTriangles.append([i,neighbour, secondNeighbour])
    triangles = np.array(listTriangles)
    vertices = np.array(listVertices)
    surfaceMesh = mesh.Mesh(np.zeros(triangles.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(triangles):
        for j in range(3):
            surfaceMesh.vectors[i][j] = vertices[f[j],:]
    surfaceMesh.save(fileName + '.stl')
    count+= 1

print("Done, " , count , " surfaces converted.")