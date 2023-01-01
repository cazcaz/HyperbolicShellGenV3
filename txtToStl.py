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
            triangleStrings = pointData[1].split("~")
            point = [float(coord) for coord in pointString]
            neighbours = [[int(index) for index in indices.split(",")] for indices in triangleStrings]
            currentSurface.append([point, neighbours])
        surfaces.append(currentSurface)
os.chdir("..")
current = os.getcwd()
os.chdir(current + "/OutputSurfaceMeshes")
for surface, fileName in zip(surfaces, fileNames):
    listVertices = []
    listTriangles = []
    point = 0
    for pointTrianglePair in surface:
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
    surfaceMesh.save(fileName + '.stl')
    count+= 1

print("Done, " , count , " surfaces converted.")