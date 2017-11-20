import sys
import math
import numpy as np
import heapq as hq
import copy
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def dot(vA, vB):
    return vA[0] * vB[0] + vA[1] * vB[1]


def ang(lineA, lineB, obtruct):
    # Get nicer vector form
    vA = [(lineA[0][0] - lineA[1][0]), (lineA[0][1] - lineA[1][1])]
    vB = [(lineB[0][0] - lineB[1][0]), (lineB[0][1] - lineB[1][1])]

    angA = math.atan2(vA[0], vA[1])
    angB = math.atan2(vB[0], vB[1])

    # Get dot prod
    dot_prod = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA) ** 0.5
    magB = dot(vB, vB) ** 0.5
    # Get cosine value
    cos_ = dot_prod / magA / magB
    # Get angle in radians and then convert to degrees
    angle = math.acos(dot_prod / magB / magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle) % 360

    if ang_deg - 180 >= 0:
        # As in if statement
        ret = 360 - ang_deg
    else:
        ret = ang_deg
    if obtruct:
        return 360 - (360 - ret), lineA[0][0], lineA[0][1]
    else:
        return 360 - ret, [lineA[0][0], lineA[0][1]]


'''
Report reflexive vertices
'''
def findReflexiveVertices(polygons):
    vertices = []
    # Your code goes here
    # You should return a list of (x,y) values as lists, i.e.
    # vertices = [[x1,y1],[x2,y2],...]
    verticesGroup = []
    for polygon in polygons:
        i = 0
        polygonArray = np.array(polygon)
        polyPath = Path(polygonArray)
        while i < len(polygon):
            if i + 1 == len(polygon):
                line_1 = [polygon[0], polygon[i]]
            else:
                line_1 = [polygon[i + 1], polygon[i]]
            if i + 2 == len(polygon):
                line_2 = [polygon[i + 1], polygon[0]]
            elif i + 1 == len(polygon):
                line_2 = [polygon[0], polygon[1]]
            else:
                line_2 = [polygon[i + 1], polygon[i + 2]]
            midpoint = [(line_2[1][0] + line_1[1][0]) / 2.0, (line_2[1][1] + line_1[1][1]) / 2.0]
            i += 1
            verticesGroup.append(ang(line_1, line_2, not polyPath.contains_point(midpoint)))

    for vertex in list(filter(lambda x: x[0] > 180, verticesGroup)):
        vertices.append(vertex[1])
    return vertices


def getExtrapoledLine(p1, p2):
    'Creates a line extrapoled in p1->p2 direction'
    xDiff = p2[0] - p1[0]
    yDiff = p2[1] - p1[1]
    distance = math.sqrt(xDiff ** 2 + yDiff ** 2)
    extendLength = 0.05
    p1_plus = [p1[0] + (p2[0] - p1[0]) / distance * extendLength, p1[1] + (p2[1] - p1[1]) / distance * extendLength]
    p1_minus = [p1[0] - (p2[0] - p1[0]) / distance * extendLength, p1[1] - (p2[1] - p1[1]) / distance * extendLength]
    p2_plus = [p2[0] + (p2[0] - p1[0]) / distance * extendLength, p2[1] + (p2[1] - p1[1]) / distance * extendLength]
    p2_minus = [p2[0] - (p2[0] - p1[0]) / distance * extendLength, p2[1] - (p2[1] - p1[1]) / distance * extendLength]
    # print plus,minus
    return p1_plus,p1_minus,p2_plus,p2_minus
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
    i=0
    while i < reflexVertices.__len__():
        vertexMap[i+1]=reflexVertices[i]
        i+=1
    edgeGroup=[]
    for polygon in polygons:
        i=0
        while i < len(polygon):
            if i == len(polygon)-1:
                edgeGroup.append((polygon[i],polygon[0]))
            else:
                edgeGroup.append((polygon[i],polygon[i+1]))
            i+=1

    obstaclesPath=[]
    for polygon in polygons:
        polygonArray = np.array(polygon)
        polyPath = Path(polygonArray)
        obstaclesPath.append(polyPath)
    i = 0
    while i < reflexVertices.__len__():
        V_A = reflexVertices[i]
        for polygon in polygons:
            try:
                k = polygon.index(V_A)
                V_A_OBIndex = polygons.index(polygon)
                break
            except ValueError:
                continue
        j = 0
        roadmapArray=[]
        while j < reflexVertices.__len__():
            if j != i:  # not connect with self
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
                F1=not obstaclesPath[V_A_OBIndex].contains_point(a)
                F2=not obstaclesPath[V_A_OBIndex].contains_point(b)
                F3=not obstaclesPath[V_B_OBIndex].contains_point(c)
                F4=not obstaclesPath[V_B_OBIndex].contains_point(d)
                F5=V_A_OBIndex==V_B_OBIndex
                F6=abs(polygons[V_A_OBIndex].index(V_A)-polygons[V_B_OBIndex].index(V_B))==1
                if F1 and F2 and F3 and F4 or F5 and F6:
                    F7 = False
                    if F1 and F2 and F3 and F4:
                                for edge in edgeGroup:
                                    B2 = edge[0] == V_A or edge[1] == V_B or edge[1] == V_A or edge[0] == V_B
                                    if not B2:
                                        B1 = intersect(V_A,V_B,edge[0],edge[1])
                                        F7= F7 or B1
                    if not F7:

                        for key in vertexMap:
                            if vertexMap[key] == V_B:
                                tmpkey = key
                        distance = round(math.sqrt((V_A[0]-V_B[0])**2+(V_A[1]-V_B[1])**2),2)
                        roadmapArray.append([tmpkey,distance])
                        # print "found line",[V_A,V_B],V_A_OBIndex,V_B_OBIndex,roadmapArray
            j+=1
        adjacencyListMap[i+1]=roadmapArray
        i+=1


    return vertexMap, adjacencyListMap


'''
Perform uniform cost search
'''
class Node:
    """
    Represents a node in the graph

    Attr:
        label: label of this node
        parent: previously visited Node before reaching current one, None by default
        f: total function cost
        g: cost to reach this node from the start
    """

    def __init__(self, label):
        """
        By default, set f = 5,000
        """
        self.label = label
        self.parent = None
        self.g = 5000
        self.f = self.g

    def __eq__(self, other):
        """
        Two cells are equivalent if their labels are equivalent
        """
        if not isinstance(other, Node):
            return False

        if self.label == other.label:
            return True
        return False

    def __str__(self):
        """
        Prints out Node in the format (label, parent's label, f)
        """
        parent_str = 'None'
        if self.parent is not None:
            parent = self.parent.label

        return "({0}, parent={1}, f={2})".format(self.label, parent_str, self.f)


def retrieve_path(start, goal, nodes_dict):
    """
    Find the path leading from start to goal by working backwards from the goal

    Parameters:
    start: label for the start vertex
    goal: label for goal vertex
    nodes_dict: dictionary with labels as keys and the corresponding vertex's node as values

    Returns:
    1D array of labels to follow from start to goal
    """
    curr_node = nodes_dict[goal]
    path = [curr_node.label]  # Start at goal

    while curr_node.label != start:
        parent = curr_node.parent
        path.append(parent.label)
        curr_node = parent

    path.reverse()  # Reverse path so it starts at start and ends at goal
    return path

def update_vertex(s, neighbor_node, cost, fringe):
    """
    Update values for a neighbor based on s

    Parameters:
    s = current Node
    neighbor_node = Node of vertex next to s
    cost = cost to move from s to neighbor
    fringe = binary heap representing fringe
    nodes_dict: dictionary with labels as keys and the corresponding vertex's node as values

    Returns: None
    """
    total_cost = s.g + cost
    if total_cost < neighbor_node.g:
        neighbor_node.g = total_cost
        neighbor_node.parent = s
        if (neighbor_node.f, neighbor_node) in fringe:
            fringe.remove((neighbor_node.f, neighbor_node))  # Remove neighbor (reorganize base on new f)

        neighbor_node.f = neighbor_node.g
        hq.heappush(fringe, (neighbor_node.f, neighbor_node))  # Insert neighbor back into fringe


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

    # Create dictionary of {labels:nodes}
    labels = adjListMap.keys()[:]
    nodes = [Node(label) for label in labels]  # Create a node for every key

    nodes_dict = {}  #
    for node in nodes:
        nodes_dict[node.label] = node

    # Run search
    start_node = nodes_dict[start]
    start_node.g = 0
    start_node.f = start_node.g
    start_node.parent = start
    fringe = []
    hq.heappush(fringe, (start_node.f, start_node))  # Insert start to fringe, need to use a 2-tuple so the heapq orders based on f-value
    closed = []  # closed := empty set

    while len(fringe) != 0:  # Checking that fringe is nonempty
        (f, s) = hq.heappop(fringe)
        if s.label == goal:
            path = retrieve_path(start, goal, nodes_dict)  # Get path from start to goal
            pathLength = nodes_dict[goal].f
            return path, pathLength
        closed.append(s.label)

        # Get neighbors and costs to move to that neighbor
        edges = adjListMap[s.label]
        neighbors = [edge[0] for edge in edges]  # Labels for neighbors of s
        edge_costs = [edge[1] for edge in edges]  # Corresponding costs to move to the neighbor

        for i in range(len(neighbors)):
            neighbor = neighbors[i]
            neighbor_node = nodes_dict[neighbor]
            edge_cost = edge_costs[i]

            if neighbor not in closed:
                update_vertex(s, neighbor_node, edge_cost, fringe)

    path = None
    pathLength = -1
    return path, pathLength

'''
Augment roadmap to include start and goal
'''
def isVisible(neighbor, point, polygons):
    """
    Check if neighbor is visible from point based on the given polygons

    :param neighbor: see if this neighbor is visible from point
    :param point: point to check visibility from
    :param polygons: list of arrays which represent coordinates in clockwise direction forming polygons
    :return: True if neighbor is visible from point, False otherwise
    """
    # Find path that gets very close to, but does not touch, the neighbor point
    # Find line between the 2 points
    x1, y1 = point
    x2, y2 = neighbor
    slope = (y2 - y1 * 1.0)/(x2 - x1 * 1.0)
    intercept = y2 - slope * x2

    # Find point right before the neighbor
    closeness = 0.999999 # Percentage of how far (with respect to distance to neighbor) new point is from the given point
    close_point_x = x1 + (x2 - x1) * closeness
    close_point_y = slope * close_point_x + intercept
    close_point = (close_point_x, close_point_y)

    test_path = Path([close_point, point], [Path.MOVETO, Path.LINETO]) # See if this segment intersects boundaries of any polygon

    for polygon in polygons:
        closed_poly = copy.deepcopy(polygon)
        closed_poly.append(closed_poly[0]) # To loop back to the beginning vertex and close the path
        num_vertices = len(closed_poly)

        # Create code that draws a line around the polygon
        codes = [Path.MOVETO] # Move to initial point
        for i in range(num_vertices-2):
            codes.append(Path.LINETO) # Every code between is LINETO (draw line from prev point to current point)
        codes.append(Path.CLOSEPOLY) # Last code closes the polygon

        polygon_path = Path(closed_poly, codes)
        if test_path.intersects_path(polygon_path):
            return False
    return True

def addToMaps(point, label, polygons, vertexMap, adjListMap):
    """
    Add a given point to the vertexMap and adjListMap

    Parameters:
    :param point: point to add to the adjListMap
    :param label: label for the given point
    :param polygons: list of list of clockwise coordinates or a polygon
    :param vertexMap: map of vertex labels to their coordinates
    :param adjListMap: adjacency list for the graph

    :return: updated adjListMap
    """
    newAdjListMap = copy.deepcopy(adjListMap)
    point_adj_list = []
    x1 = point[0]
    y1 = point[1]

    # Add edges between given point and other visible points
    for vertex_label in newAdjListMap.keys():
        vertex = vertexMap[vertex_label]
        if isVisible(vertex, point, polygons):
            x2 = vertex[0]
            y2 = vertex[1]
            distance = math.sqrt(math.pow(x2 - x1, 2) + math.pow(y2 - y1, 2))
            distance = round(distance, 2)
            newAdjListMap[vertex_label].append([label, distance])
            point_adj_list.append([vertex_label, distance])

    newAdjListMap[label] = point_adj_list

    # Add point to vertex map
    newVertexMap = copy.deepcopy(vertexMap)
    newVertexMap[label] = point
    return newVertexMap, newAdjListMap

def updateRoadmap(polygons, vertexMap, adjListMap, x1, y1, x2, y2):
    updatedALMap = dict()
    startLabel = 0
    goalLabel = -1

    # Your code goes here. Note that for convenience, we
    # let start and goal have vertex labels 0 and -1,
    # respectively. Make sure you use these as your labels
    # for the start and goal vertices in the shortest path
    # roadmap. Note that what you do here is similar to
    # when you construct the roadmap.

    vertex_map_copy = copy.deepcopy(vertexMap)

    start = (x1, y1)
    goal = (x2, y2)
    vertex_map_copy, updatedALMap = addToMaps(start, startLabel, polygons, vertex_map_copy, adjListMap)
    vertex_map_copy, updatedALMap = addToMaps(goal, goalLabel, polygons, vertex_map_copy, updatedALMap)

    return startLabel, goalLabel, updatedALMap

'''
Visualize the roadmap (including start and goal) in green, path computed in red
'''
# Set up matplotlib to create a plot with an empty square (from visualize.py)
def setupPlot():
    fig = plt.figure(num=None, figsize=(5, 5), dpi=120, facecolor='w', edgecolor='k')
    plt.autoscale(False)
    plt.axis('off')
    ax = fig.add_subplot(1,1,1)
    ax.set_axis_off()
    ax.add_patch(patches.Rectangle(
        (0,0),   # (x,y)
        1,          # width
        1,          # height
        fill=False
        ))
    return fig, ax

# Make a patch for a single poly (from visualize.py)
def createPolygonPatch(polygon):
    verts = []
    codes= []
    for v in range(0, len(polygon)):
        xy = polygon[v]
        verts.append((xy[0]/10., xy[1]/10.))
        if v == 0:
            codes.append(Path.MOVETO)
        else:
            codes.append(Path.LINETO)
    verts.append(verts[0])
    codes.append(Path.CLOSEPOLY)
    test_path = Path(verts, codes)
    patch = patches.PathPatch(test_path, facecolor='gray', lw=1)

    return patch

# Make a patch for the robot (from visualize.py)
def createPolygonPatchForRobot(polygon):
    verts = []
    codes= []
    for v in range(0, len(polygon)):
        xy = polygon[v]
        verts.append((xy[0]/10., xy[1]/10.))
        if v == 0:
            codes.append(Path.MOVETO)
        else:
            codes.append(Path.LINETO)
    verts.append(verts[0])
    codes.append(Path.CLOSEPOLY)
    test_path = Path(verts, codes)
    patch = patches.PathPatch(test_path, facecolor='gray', lw=1)

    return patch

def visualize(file_name, polygons, vertexMap, adjListMap, start, x1, y1, goal, x2, y2, path):
    """
    Plot the given roadmap, obstacles and path along it

    Parameters:
    file_name = name of the file with obstacle information
    polygons = list of polygons (their vertices in clockwise order) that act as obstacles
    vertexMap = map from vertex labels to their coordinates
    adjListMap = adjacency list mapping vertices to their outgoing edges
    start = label for start vertex
    x1 = x-coordinate of start
    y1 = y-coordiante of start
    goal = label for goal vertex
    x2 = x-coordinate of goal
    y2 = y-coordinate of goal
    path = labels of vertices that lead from start to goal

    Returns: None
    """
    # From visualize.py
    fig, ax = setupPlot()
    for p in range(0, len(polygons)):
        patch = createPolygonPatch(polygons[p])
        ax.add_patch(patch)
    # From visualize.py

    # Create vertex map that includes start and goal vertices
    updatedVertexMap = copy.deepcopy(vertexMap)
    updatedVertexMap[start] = (x1, y1)
    updatedVertexMap[goal] = (x2, y2)

    # Normalize all vertex coordinates (because polygons are outputted in 1x1 box, not 10x10 box)
    for label in updatedVertexMap.keys():
        vertex_pt = updatedVertexMap[label]
        updatedVertexMap[label] = [vertex_pt[0]/10.0, vertex_pt[1]/10.0]

    # Add roadmap to plot
    roadmap_x = []
    roadmap_y = []
    for vertex in adjListMap.keys():
        vertex_pt = updatedVertexMap[vertex]
        roadmap_x = []
        roadmap_y = []

        # Find all neighbors labels
        neighbors = [x[0] for x in adjListMap[vertex]]

        # Draw green line segment from inital vertex to its neighors
        for neighbor in neighbors:
            neighbor_pt = updatedVertexMap[neighbor]
            roadmap_x = [vertex_pt[0]/10.0, neighbor_pt[0]/10.0]
            roadmap_y = [vertex_pt[1]/10.0, neighbor_pt[1]/10.0]
            plt.plot(roadmap_x, roadmap_y, 'g')

    # Add path to plot
    path_x = []
    path_y = []

    for label in path:
        pt = updatedVertexMap[label]
        path_x.append(pt[0])
        path_y.append(pt[1])
    plt.plot(path_x, path_y, 'r')

    # Plot data
    title = "Shortest Roadmap ({})".format(file_name)
    plt.title(title)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":

    # Retrive file name for input data
    if (len(sys.argv) < 6):
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
    start, goal, updatedALMap = updateRoadmap(polygons,vertexMap,adjListMap, x1, y1, x2, y2)
    print "Updated roadmap:"
    print str(updatedALMap)
    print ""

    # Search for a solution
    path, length = uniformCostSearch(updatedALMap, start, goal)
    print "Final path:"
    print str(path)
    print "Final path length:" + str(length)

    # Extra visualization elements goes here
    visualize(filename, polygons, vertexMap, updatedALMap, start, x1, y1, goal, x2, y2, path)

