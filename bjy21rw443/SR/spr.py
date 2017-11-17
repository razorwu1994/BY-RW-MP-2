import sys
import math
import numpy as np
import matplotlib.path as path

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]
def ang(lineA, lineB, obtruct):
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]

    angA = math.atan2(vA[0],vA[1])
    angB = math.atan2(vB[0],vB[1])

    # Get dot prod
    dot_prod = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = dot_prod/magA/magB
    # Get angle in radians and then convert to degrees
    angle = math.acos(dot_prod/magB/magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360

    if ang_deg-180>=0:
        # As in if statement
        ret = 360 - ang_deg
    else:
        ret = ang_deg
    if obtruct:
        return 360-(360-ret),lineA[0][0],lineA[0][1]
    else:
        return 360-ret,[lineA[0][0],lineA[0][1]]

'''
Report reflexive vertices
'''
def findReflexiveVertices(polygons):
    vertices=[]
    # Your code goes here
    # You should return a list of (x,y) values as lists, i.e.
    # vertices = [[x1,y1],[x2,y2],...]
    verticesGroup=[]
    for polygon in polygons:
        i=0
        polygonArray = np.array(polygon)
        polyPath = path.Path(polygonArray)
        while i < len(polygon):
            if i+1==len(polygon):
                line_1=[polygon[0],polygon[i]]
            else:
                line_1=[polygon[i+1],polygon[i]]
            if i+2==len(polygon):
                line_2=[polygon[i+1],polygon[0]]
            elif i+1==len(polygon):
                line_2=[polygon[0],polygon[1]]
            else:
                line_2=[polygon[i+1],polygon[i+2]]
            midpoint = [(line_2[1][0]+line_1[1][0])/2.0,(line_2[1][1]+line_1[1][1])/2.0]
            i+=1
            verticesGroup.append(ang(line_1,line_2,not polyPath.contains_point(midpoint)))

    for vertex in list(filter(lambda x: x[0] > 180, verticesGroup)):
        vertices.append(vertex[1])
    return vertices

def getExtrapoledLine(p1,p2):
    'Creates a line extrapoled in p1->p2 direction'
    xDiff = p2[0]-p1[0]
    yDiff = p2[1]-p1[1]
    distance = math.sqrt(xDiff**2+yDiff**2)
    extendLength=0.05
    p1_plus = [p1[0]+(p2[0]-p1[0])/distance*extendLength,p1[1]+(p2[1]-p1[1])/distance*extendLength]
    p1_minus = [p1[0]-(p2[0]-p1[0])/distance*extendLength,p1[1]-(p2[1]-p1[1])/distance*extendLength]
    p2_plus = [p2[0]+(p2[0]-p1[0])/distance*extendLength,p2[1]+(p2[1]-p1[1])/distance*extendLength]
    p2_minus = [p2[0]-(p2[0]-p1[0])/distance*extendLength,p2[1]-(p2[1]-p1[1])/distance*extendLength]
    # print plus,minus
    return p1_plus,p1_minus,p2_plus,p2_minus
    #print p1,path.Path(np.array([p1,p2])).contains_point(p1),p2,path.Path(np.array([p1,p2])).contains_point(p2)
    # return a,b
'''
Compute the roadmap graph
'''
def computeSPRoadmap(polygons, reflexVertices):
    vertexMap = dict()
    adjacencyListMap = dict()

    # Your code goes here
    # You should check for each pair of vertices whether the
    # edge between them should belong to the shortest path
    # roadmap.
    #
    # Your vertexMap should look like
    # {1: [5.2,6.7], 2: [9.2,2.3], ... }
    #
    # and your adjacencyListMap should look like
    # {1: [[2, 5.95], [3, 4.72]], 2: [[1, 5.95], [5,3.52]], ... }
    #
    # The vertex labels used here should start from 1
    obstaclesPath=[]
    for polygon in polygons:
        polygonArray = np.array(polygon)
        polyPath = path.Path(polygonArray)
        obstaclesPath.append(polyPath)
    i=0
    while i < reflexVertices.__len__():
        V_A = reflexVertices[i]
        for polygon in polygons:
            try:
                k = polygon.index(V_A)
                V_A_OBIndex = polygons.index(polygon)
                break
            except ValueError:
                continue
        j=0
        while j < reflexVertices.__len__():
            if j!=i:#not connect with self
                V_B = reflexVertices[j]
                V_B_OBIndex = V_A_OBIndex
                for polygon in polygons:
                    try:
                        k = polygon.index(V_B)
                        V_B_OBIndex = polygons.index(polygon)
                        break
                    except ValueError:
                        continue
                a,b,c,d = getExtrapoledLine(V_A,V_B)
                # midpoint = [(V_A[0]+V_B[0])/2.0,(V_A[1]+V_B[1])/2.0]
                # print not obstaclesPath[V_B_OBIndex].contains_point(a),not obstaclesPath[V_B_OBIndex].contains_point(b)
                F1=not obstaclesPath[V_A_OBIndex].contains_point(a)
                F2=not obstaclesPath[V_A_OBIndex].contains_point(b)
                F3=not obstaclesPath[V_B_OBIndex].contains_point(c)
                F4=not obstaclesPath[V_B_OBIndex].contains_point(d)
                F5=V_A_OBIndex==V_B_OBIndex
                F6=abs(polygons[V_A_OBIndex].index(V_A)-polygons[V_B_OBIndex].index(V_B))==1
                if F1 and F2 and F3 and F4 or F5 and F6:
                    if F5 and F6:
                        print "same poly"
                    print "found line",[V_A,V_B],V_A_OBIndex,V_B_OBIndex
            j+=1
        i+=1


    return vertexMap, adjacencyListMap

'''
Perform uniform cost search
'''
def uniformCostSearch(adjListMap, start, goal):
    path = []
    pathLength = 0

    # Your code goes here. As the result, the function should
    # return a list of vertex labels, e.g.
    #
    # path = [23, 15, 9, ..., 37]
    #
    # in which 23 would be the label for the start and 37 the
    # label for the goal.

    return path, pathLength

'''
Agument roadmap to include start and goal
'''
def updateRoadmap(adjListMap, x1, y1, x2, y2):
    updatedALMap = dict()
    startLabel = 0
    goalLabel = -1

    # Your code goes here. Note that for convenience, we
    # let start and goal have vertex labels 0 and -1,
    # respectively. Make sure you use these as your labels
    # for the start and goal vertices in the shortest path
    # roadmap. Note that what you do here is similar to
    # when you construct the roadmap.

    return startLabel, goalLabel, updatedALMap

if __name__ == "__main__":

    # Retrive file name for input data
    if(len(sys.argv) < 6):
        print "Five arguments required: python spr.py [env-file] [x1] [y1] [x2] [y2]"
        exit()

    filename = sys.argv[1]
    x1 = float(sys.argv[2])
    y1 = float(sys.argv[3])
    x2 = float(sys.argv[4])
    y2 = float(sys.argv[5])

    # Read data and parse polygons
    lines = [line.rstrip('\n') for line in open(filename)]
    polygons = []
    for line in range(0, len(lines)):
        xys = lines[line].split(';')
        polygon = []
        for p in range(0, len(xys)):
            polygon.append(map(float, xys[p].split(',')))
        polygons.append(polygon)

    # Print out the data
    print "Pologonal obstacles:"
    for p in range(0, len(polygons)):
        print str(polygons[p])
    print ""

    # Compute reflex vertices
    reflexVertices = findReflexiveVertices(polygons)
    print "Reflexive vertices:"
    print str(reflexVertices)
    print ""

    # Compute the roadmap
    vertexMap, adjListMap = computeSPRoadmap(polygons, reflexVertices)
    print "Vertex map:"
    print str(vertexMap)
    print ""
    print "Base roadmap:"
    print str(adjListMap)
    print ""

    # Update roadmap
    start, goal, updatedALMap = updateRoadmap(adjListMap, x1, y1, x2, y2)
    print "Updated roadmap:"
    print str(updatedALMap)
    print ""

    # Search for a solution
    path, length = uniformCostSearch(updatedALMap, start, goal)
    print "Final path:"
    print str(path)
    print "Final path length:" + str(length)


    # Extra visualization elements goes here