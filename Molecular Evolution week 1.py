#Week 1
#Exercise 1
import sys
import re #proides regular expressions which help extract numbers from strings
import queue

WeightedEdge = tuple[int, int, int] #each edge contains: starting node, ending node, edge weight)
DistMatrix = list[list[int]] #distance matrix is a 2D list containing integer distances


def distances_between_leaves(n: int, edge_list: list[WeightedEdge]) -> DistMatrix:
    # 1. Parse edges safely (handles both tuples and string formatted edges)
    parsed_edges: list[tuple[int, int, int]] = [] #create an empty list to store all edges in a standardized format
    max_node_id = -1 #maximum node ID is set to -1, to keep track of larhest node number encountered

    for item in edge_list:
        if isinstance(item, str): #check if current edge is represented as a string
            match = re.search(r"(\d+)(?:->|\s+)(\d+)(?::|\s+)(\d+)", item.strip()) #use regular exprpession to extract 3 numbers
            if match:
                u, v, w = map(int, match.groups()) #extract 3 captured values and convert them to integers
                parsed_edges.append((u, v, w))
                max_node_id = max(max_node_id, u, v)
        elif isinstance(item, (list, tuple)) and len(item) == 3: #otherwise check if edge is represented as a list or a tuple and that it contains 3 elements
            u, v, w = int(item[0]), int(item[1]), int(item[2]) #extract the starting node, ending node and edge weight, convert all 3 values to integers
            parsed_edges.append((u, v, w)) #add edge to the standardized edge list as a tuple
            max_node_id = max(max_node_id, u, v) #update maximum node ID using the 2 node numbers, if edge weight is not a node ID it is excluded

    # Base cases
    if max_node_id == -1:
        return [[0]] if n == 1 else [] #if there is only one leaf, distance itself is 0, otherwise return an empty list

    # 2. Allocate adjacency list and degree array dynamically based on the max node ID found
    total_nodes = max(n, max_node_id + 1) #calculate total number of nodes needed
    #max_id + 1 is used because python indexing starts at 0
    adj: list[list[tuple[int, int]]] = [[] for _ in range(total_nodes)] #create an empty adjacency list for all noes, each node will store its neighbors and the corresponding edge weight
    degrees: list[int] = [0] * total_nodes #create a list to store a degree (number of edges) of every node, initially every node has degree 0

    for u, v, w in parsed_edges:#all standardized edges
        adj[u].append((v, w)) #add node v as a neighbor of node u with edge w
        adj[v].append((u, w)) #add node u as a neighbor of node v with edge w
        degrees[u] += 1
        degrees[v] += 1

    # 3. Identify leaves (nodes with degree 1) sorted by node ID
    leaves: list[int] = [i for i in range(total_nodes) if degrees[i] == 1] #check every node and select those with degree 1
    leaves.sort()

    leaf_count = len(leaves)
    if leaf_count == 0:
        return []

    dist_matrix: DistMatrix = [[0] * leaf_count for _ in range(leaf_count)] #distance matrix with number of leaves x number of leaves, initially all distances set to 0
    leaf_index_map: dict[int, int] = {leaf: idx for idx, leaf in enumerate(leaves)} #a dictionary that maps each leaf;s original node ID to its corresponding index in distance matrix

    # 4. Traversal from each leaf to measure distances to all other leaves
    for start_idx, start_node in enumerate(leaves): #start_idx represents the leaf's psoition in distance matrix, start_node represents its actual node ID in tree
        stack = [(start_node, -1, 0)] #each stack element contains 3 values, we start at start_node, with no parent (-1) and distance = 0
        while stack:
            curr, parent, dist = stack.pop() #remove the last element from stack, extract current node, its parent and accumulated distance

            if curr in leaf_index_map: #check if current node is a leaf
                target_idx = leaf_index_map[curr] #find position of current leaf in distance matrix
                dist_matrix[start_idx][target_idx] = dist #store the calculated distance between starting leaf and current leaf in appropriate matrix position

            for neighbor, weight in adj[curr]: #explore all neighbors of current node, each neighbor is repesented by (neighbor node, node weight)
                if neighbor != parent: #avoid going back to parent node because tree is unidirected 
                    stack.append((neighbor, curr, dist + weight)) #add neighbor node to stack,neighbor becomes next current node

    return dist_matrix
  
#Exercise 2
from math import inf #used to initialize minimum limb length to a very large value

DistMatrix = list[list[float]]

def limb_length(j: int, distance_matrix: DistMatrix) -> float: #calculate limb length of leaf such as j is the index of leaf whose limb length we want to calculate
    """
    Calculates the limb length of leaf j from an n x n distance matrix.
    
    Handles base cases (n <= 2) to prevent returning inf.
    """
    n = len(distance_matrix) #gets number of rows in distance matrix (and number of columns because  matrix is square)
    
    # Base cases for small matrices
    if n <= 1: #if we have 0 or 1 leaves teh return 0
        return 0.0
    if n == 2:
        val = distance_matrix[0][1] #if there are exactly 2 leaves return distance between them
        return int(val) if isinstance(val, float) and val.is_integer() else val #check if value is a float (5.0) representing a whole number (5), if so convert it to an integer

    min_limb = inf
    
    # Standard formula for n >= 3: min (D[i,j] + D[j,k] - D[i,k]) / 2
    for i in range(n):
        if i == j:
            continue #skip
        for k in range(i + 1, n): #to avoid checking duplicate pairs
            if k == j:
                continue #skip
            
            current_limb = (distance_matrix[i][j] + distance_matrix[j][k] - distance_matrix[i][k]) / 2
            if current_limb < min_limb: #compare newly calculated candidate with smalles limb length so far
                min_limb = current_limb #if candidate is smaller update min_limb
                
    return int(min_limb) if isinstance(min_limb, float) and min_limb.is_integer() else min_limb
  
#Exercise 3
import collections

DistMatrix = list[list[int]] #matrix contains distance between pairs of leaves
WeightedEdge = tuple[int, int, float] #each edge contains a starting node, an ending noe and edge weight

def limb_length(j: int, distance_matrix: DistMatrix, current_n: int) -> float: #distance matrix contains distances between all leaves
    """Calculates limb length of leaf j restricted to current_n nodes.""" #current_n indicates how many leaves are currently being conisdered
    if current_n <= 1:
        return 0.0
    if current_n == 2:
        return distance_matrix[0][1]

    min_limb = float('inf')
    for i in range(current_n):
        if i == j:
            continue
        for k in range(i + 1, current_n):
            if k == j:
                continue
            limb = (distance_matrix[i][j] + distance_matrix[j][k] - distance_matrix[i][k]) / 2
            if limb < min_limb:
                min_limb = limb
    return min_limb


def additive_phylogeny(distance_matrix: DistMatrix) -> list[WeightedEdge]: #receives a distance matrix and recinstructs a weighted tree whose leaf to leaf distances match the given matrix
    n = len(distance_matrix)
    next_node = n

    def tree_to_adj(adj: dict[int, dict[int, float]]) -> list[WeightedEdge]: #a function to convert the adjacency dictionary into a list of weighted edges
        edges = []
        for u in adj:
            for v, w in adj[u].items(): #iterate through all neighbors of node u and their edge weights
                edges.append((u, v, int(w) if float(w).is_integer() else w)) #add the edge (u, v, w) to output list, if weight is a whole-number floar convert it to an integer
        return edges

    def add_edge(adj: dict, u: int, v: int, w: float): #adds an undirected weighted edge between nodes u and v
        adj[u][v] = w #add v as a neighbor of u with weight w
        adj[v][u] = w #add u as a neighbor of v with same weight w

    def remove_edge(adj: dict, u: int, v: int): #removes undirected edge between u and v
        adj[u].pop(v, None) #removes v from neighbors of u (None is defaukt prevents error if edge does not exist)
        adj[v].pop(u, None) #removes u from neighbors of v

    def find_path(adj: dict, start: int, target: int) -> list[int]: #finds path between 2 nodes
        queue = collections.deque([[start]]) #queue containing starting path, each queue element represents a complete path
        #deque removes elements from the front
        visited = {start} #creates a set to keep track of visited nodes to avoid revisiting them
        while queue:
            path = queue.popleft() #remoe first path from queue, first in first out principle
            curr = path[-1] #gets the last node, which we are currently exploring
            if curr == target:
                return path
            for neighbor in adj[curr]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(path + [neighbor]) #add a new path to queue which consists of previous path  newly discovered neighbor
        return []

    def solve(matrix: DistMatrix, current_n: int) -> dict[int, dict[int, float]]: #reconstructs the tree
        nonlocal next_node #we want to use and modify the next_node variable defined in outer additive_phylogeny function to ensure that every new internal node gets a unique ID

        # Base case: 2 leaves
        if current_n == 2:
            adj = collections.defaultdict(dict) #creates an adjacency dictionary with an empty neighbor dictionary whenever a new node is accessed
            add_edge(adj, 0, 1, matrix[0][1]) #connect leaf 1 to leaf 2, edge weight is the distance between these 2 leaves
            return adj

        leaf = current_n - 1
        
        # 1. Calculate limb length using current un-mutated submatrix state
        limb = limb_length(leaf, matrix, current_n)

        # 2. Find attachments leaves i and k in original submatrix, will store teh 2 leaves between which the selected leaf must be attached
        #-1 indicates that no valid pair has been found yet
        i_found, k_found = -1, -1
        x_dist = 0 #x-dist represents the distance from leaf i_found to attachement point of removed leaf
        for i in range(current_n - 1): #examine all leaves except leaf currently being removed
            for k in range(current_n - 1): #examine all possible second leaves except the removed leaf
                if i != k and matrix[i][k] == (matrix[i][leaf] - limb) + (matrix[leaf][k] - limb): #if thsi condirion is satisfied the attachement point of leaf j lies on path between leaves i and k
                    i_found, k_found = i, k #storing first and second attachement leaves
                    x_dist = matrix[i][leaf] - limb #distance from leaf i to attachement point = distance from i to j - j's limb length
                    break
            if i_found != -1:
                break

        # 3. Create a trimmed copy for the recursive call
        trimmed_matrix = [row[:current_n - 1] for row in matrix[:current_n - 1]] #a smaller disnatnce matrix containing only first current_n -1 leaves, no longer contains the removed leaf
        
        # 4. Recursively build tree for first current_n - 1 leaves, each call removes one leaf until only 2 leaves remain
        adj = solve(trimmed_matrix, current_n - 1)

        # 5. Attach leaf back into the tree
        path = find_path(adj, i_found, k_found) #find unique path between the 2 attachement leaves in reconstructed tree
        dist_accumulated = 0.0 #initialize the accumulated distance from first arttachement leaf, this variable tracks how far we have traveled along the path

        for idx in range(len(path) - 1):# traverses th epath one edge at a time, len(path) - 1 is used because each edge connects 2 consecutive nodes
            u, v = path[idx], path[idx + 1] #gets the 2 endpoints of current edge, u is the current node and v is the next node on the path
            edge_weight = adj[u][v] #retried the weight of edge connecting u and v

            if dist_accumulated + edge_weight < x_dist: #if attachement point lies beyond current edge (= we have not reached x_dist)
                dist_accumulated += edge_weight #add current edge weight to accumulated distance
            elif dist_accumulated + edge_weight == x_dist: #if attacahement node is exactly at node v, the accumulated distance after traversing the edge = x_dist
                add_edge(adj, leaf, v, limb) #attached the removed leaf directly to node v, the edge weight is the limb lenght of removed leaf
                break
            else: #if attachement point lies somewhere inside current edge, we must plit edge and create an internal node
                remove_edge(adj, u, v) #remove original edge between u and v

                new_node = next_node #assign next available node ID to new internal node
                next_node += 1 #increment next_node so the next internal node gets a unique ID

                dist_to_new = x_dist - dist_accumulated #calculate the distance from u to internal node, this is the remaining distance needed to reach x_dist
                dist_from_new = edge_weight - dist_to_new #calculate the remining distance from new internal node to original edge weight is divided into 2 parts

                add_edge(adj, u, new_node, dist_to_new) #connect u to nrely created internal node, the edge weight is distane from u to new node
                add_edge(adj, new_node, v, dist_from_new)b#edge weight is the remaining portion of original edge
                add_edge(adj, new_node, leaf, limb) #attach removed leaf to new internal node, the connecting edge has the calculated limb length
                break

        return adj

    return tree_to_adj(solve(distance_matrix, n))

