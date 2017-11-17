import sys
import math
import numpy as np
import matplotlib.path as path
import heapq as hq

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
class Node:
    """
    Represents a node in the graph

    Attr:
        pos: coordinates (x, y) of this node in  C-space
        parent: previously visited Cell before reaching current one, None by default
        f: function value, equivalent to g+h
    """

    def __init__(self, pos):
        """
        By default, set f = 20,000
        """
        self.pos = pos
        self.parent = None
        self.f = 20000

    def __eq__(self, other):
        """
        Compare two cells based on their f (priority) values
        """
        if not isinstance(other, Cell):
            return False

        if self.pos == other.pos:
            return True
        return False

    def __str__(self):
        """
        Prints out the card with the format: <value> of <suit>
        Jokers are just printed out as 'joker'
        """
        t_type = self.convert_to_char()
        return "(({0}, {1}), {2}, f={3}, g={4}, h={5})".format(self.pos[0], self.pos[1], t_type, self.f, self.g, self.h)

def retrieve_path(start, goal, grid):
    """
    Find the path leading from start to goal by working backwards from the goal

    Parameters:
    start: (x, y) coordinates of the start position
    goal: (x, y) coordaintes of goal position
    grid: 160x120 array of Cells

    Returns:
    1D array of (x, y) coordinates to follow from start to goal
    """
    curr_cell = grid[goal[0]][goal[1]]
    path = [curr_cell.pos]  # Start at goal

    while curr_cell.pos != start:
        parent = curr_cell.parent
        path.append(parent.pos)
        curr_cell = parent

    path.reverse()  # Reverse path so it starts at start and ends at goal
    return path

def get_neighbors(cell, grid):
    """
    Find the valid neighbors for the given cell.
    Check 8-neighbors around the cell, ignore blocked cells and cells outside of the boundary.

    Parameters:
    cell = target Cell
    grid = 160x120 grid of Cells

    Returns: 1D array of Cells
    """
    # Find 8 neighboring positions
    pos = cell.pos

    top_left_pos = (pos[0] - 1, pos[1] + 1)
    top_pos = (pos[0], pos[1] + 1)
    top_right_pos = (pos[0] + 1, pos[1] + 1)
    right_pos = (pos[0] + 1, pos[1])
    bottom_right_pos = (pos[0] + 1, pos[1] - 1)
    bottom_pos = (pos[0], pos[1] - 1)
    bottom_left_pos = (pos[0] - 1, pos[1] - 1)
    left_pos = (pos[0] - 1, pos[1])

    possible_neighbors = [top_left_pos, top_pos, top_right_pos, right_pos, bottom_right_pos, bottom_pos,
                          bottom_left_pos, left_pos]

    # Filter out invalid neighbors (out of bounds or blocked cell)
    possible_neighbors = [pos for pos in possible_neighbors if
                          pos[0] >= 0 and pos[0] < 120 and pos[1] >= 0 and pos[1] < 160]

    for neighbor in possible_neighbors:
        if grid[neighbor[0]][neighbor[1]].terrain_type == BLOCKED:
            possible_neighbors.remove(neighbor)

    """ Testing
    print "Neighbors:"
    for neighbor in possible_neighbors:
        print neighbor
    print ""
    """

    valid_neighbors = [grid[pos[0]][pos[1]] for pos in possible_neighbors]
    return valid_neighbors

def get_cost(s, neighbor):
    """
    Calculate cost to move from s to its neighboring cell.
        - Unblocked cells cost 1 to traverse along an edge
        - Hard-to-traverse cells cost 2 to traverse along an edge
        - Moving from highway to highway cuts overall cost by a factor of 4
    Parameter:
    s = Cell for the furthest cell on the optimal path
    neighbor = Cell for a neighbor of s

    Returns: cost to move from s to neighbor
    """
    # Find Euclidean distance
    (x1, y1) = s.pos
    (x2, y2) = neighbor.pos
    distance = math.sqrt(math.pow(x2 - x1, 2) + math.pow(y2 - y1, 2))

    # Factor in terrain types
    temp_cells = [s, neighbor]
    temp_dists = [distance / 2, distance / 2]

    for i in range(len(temp_cells)):
        if temp_cells[i].terrain_type == ROUGH:
            temp_dists[i] *= 2  # Rough terrain costs 2x to move across

    distance = sum(temp_dists)

    # Check if highway cuts cost further (4x)
    if s.has_highway is True and neighbor.has_highway is True:
        distance /= 4

    return distance

def get_heuristic(cell, grid):
    """
    Calculate the heursitic for a cell

    Parameters:
    cell = target cell
    grid = 160x120 grid

    Returns: h value for the cell
    """
    return 0  # For UCS use 0, replace with something else for A* and weighted A*

def update_vertex(s, neighbor, fringe):
    """
    Update values for a neighbor based on s

    Parameters:
    s = a Cell
    neighbor = a Cell next to s

    Returns: None
    """
    total_cost = s.g + get_cost(s, neighbor)
    if total_cost < neighbor.g:
        neighbor.g = total_cost
        neighbor.parent = s
        if (neighbor.f, neighbor) in fringe:
            fringe.remove((neighbor.f, neighbor))  # Possible optimization opportunity?

        neighbor.f = neighbor.g + neighbor.h  # Update neighbor's f-value
        hq.heappush(fringe, (neighbor.f, neighbor))  # Insert neighbor back into fringe

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

    # Run search
    start_cell = grid[start[0]][start[1]]
    start_cell.g = 0
    start_cell.f = start_cell.g + start_cell.h
    start_cell.parent = start
    fringe = []
    hq.heappush(fringe, (start_cell.f, start_cell)) # Insert start to fringe, need to use a 2-tuple so the heapq orders based on f-value
    closed = [] # closed := empty set

    while len(fringe) != 0: # Checking that fringe is nonempty
        (f, s) = hq.heappop(fringe)
        if s.pos == goal:
            path = retrieve_path(start, goal, grid) # Get path from start to goal
            pathLength = len(path)
            return path, pathLength
        closed.append(s.pos)
        neighbors = get_neighbors(s, grid)
        for neighbor in neighbors:
            if neighbor.pos not in closed: # Possible optimization opportunity
                if (neighbor.f, neighbor) not in fringe:
                    neighbor.g = 20000 # 20,000 = infinity
                    neighbor.parent = None
                update_vertex(s, neighbor, fringe)
    path = None
    pathLength = 0
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
