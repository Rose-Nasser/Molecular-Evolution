#Week 2
#Exercise 1
import sys


def parse_input():
  """Parses the distance matrix from standard input.
    
    Handles inputs separated by tabs or spaces, with or without an explicit dimension n on the first line.
    """
raw_lines = [line.strip() for line in sys.stdin.read().strip().splitlines() if line.strip()] #removes all input from standard input
#strip() removes unnecessary spaces and newline characters at the beginning and end
#splitlines() separates the input into individual lines
if not raw_lines:
  return []

# Split on ANY consecutive whitespace (tabs or spaces)
first_line_tokens = raw_lines[0].split()

# If the first line contains a single integer, it specifies 'n'
if len(first_line_tokens) == 1:
  n = int(first_line_tokens[0]) # this is the number of rows and columns in distance matrix
matrix_lines = raw_lines[1:n + 1] #exrract the next n lines because the first line contains the dimensions not a row of the distance matrix
else:
  # Matrix starts directly on line 0
  n = len(raw_lines) #if first line contains multiple values, assume the distance matrix starts directly on first line. In this case, the number of matrix rows equals the number of input lines
matrix_lines = raw_lines[:n] #select the first n lines as the matrix rows

matrix = [] #create an empty list to store the complete distance matrix
for line in matrix_lines:
  row = [float(x) for x in line.split()] #split current line into individual values using whitespace and convert any value into float
matrix.append(row)

return matrix


def upgma(D): #unweighted pair groups method with arithmetic mean
  """Executes the UPGMA algorithm to construct an ultrametric tree.

    Args:
        D: n x n distance matrix (list of lists of floats/ints).

    Returns:
        adj: Dict mapping parent node ID -> list of (child_id, edge_weight).
    """
n = len(D) #calculate the number of leaves in original distance matrix
if n == 0:
  return {}

# Active clusters mapping: cluster_id -> {"size": int, "age": float}
clusters = {i: {"size": 1, "age": 0.0} for i in range(n)} #each cluster has two properties: size (number of original leaves contained in cluster) and age (height of cluster in phylogenetic tree)

# Distance map: (cluster_i, cluster_j) -> distance
dist = {}
for i in range(n):
  for j in range(n):
  dist[(i, j)] = float(D[i][j]) #store the distance between cluster i and cluster j. Convert the value to float to ensure consistency

adj = {}
next_node_id = n #assign the next available internal node ID. Original leaves have IDs from 0 to n-1

# Main clustering loop
while len(clusters) > 1: #continue merging clusters until only one cluster remains
  min_dist = float("inf")
closest_pair = (-1, -1) #initialize the closest pair with placeholder values, these values will be replaced when the closest pair is found.

cluster_ids = sorted(clusters.keys()) #get ids of all currently active clusters
for i in range(len(cluster_ids)): #outer loop selects first cluster
  for j in range(i + 1, len(cluster_ids)): #inner loop selects second cluster,m starting from i + 1 avoids comparing a cluster with itself and prevents checking the same pair twice
  c1, c2 = cluster_ids[i], cluster_ids[j] #retrieve ids of two clusters being compared
if dist[(c1, c2)] < min_dist: #check whether the distance between these clusters is smaller than minimum distance found so far
  min_dist = dist[(c1, c2)] #update minimum distance
closest_pair = (c1, c2) #store ids of the two closest clusters

c_i, c_j = closest_pair #retrieve the ids of the two clusters that will be merged
c_new = next_node_id #assign a new traditional node ID to represent the merged cluster
next_node_id += 1 #increment the counter so that next internal node gets a unique ID

# Calculate age of new cluster node
age_new = min_dist / 2.0 #in upgma the age of newly created cluster is half the distance between two clusters being merged
size_new = clusters[c_i]["size"] + clusters[c_j]["size"] #calculate the size of the new cluster, the size is the sum of sizes of two merged clusters
clusters[c_new] = {"size": size_new, "age": age_new} #add the new cluster to the active clusters dictionary. Store its total size and its calculated age

# Calculate directed edge lengths from parent (c_new) to children (c_i, c_j)
len_i = age_new - clusters[c_i]["age"]
len_j = age_new - clusters[c_j]["age"]

# Add directed edges from parent to children
if c_new not in adj: #checks  whether the new internal node already exists in adjacency dictionary
  adj[c_new] = [] #if it does not exist create an empty list for its children
adj[c_new].append((c_i, len_i))
adj[c_new].append((c_j, len_j))

# Compute new distances to all remaining active clusters
for c in clusters: #iterate through all clusters currently stored in the dictionary, including the two clusters that were just merged and the new cluster
  if c != c_i and c != c_j and c != c_new: #skip the two merged clusters and the newly created cluster. We will calculate distances between the new cluster and the other emaining active clusters
  d_i = dist[(c_i, c)] #retrieve the distance between the first emerged cluster and the remaining cluster c
d_j = dist[(c_j, c)] #retrieve the distance between the second merged cluster and the remaining clustr c
s_i = clusters[c_i]["size"] #retrieve the number of original leaves in the first cluster
s_j = clusters[c_j]["size"] #retrieve the number of original leaves in the second cluster

# Weighted arithmetic average
d_new = (d_i * s_i + d_j * s_j) / (s_i + s_j)
dist[(c_new, c)] = d_new #store the distance from new cluster to cluster c
dist[(c, c_new)] = d_new #store the reverse distance as well because distance matrix is symmetric

# Remove first and second merged clusters from active set
del clusters[c_i]
del clusters[c_j]

return adj


def format_weight(weight):
  """Formats edge weights to match grader output specs:
    - If weight is an integer (e.g. 5.0), prints '5'
    - Otherwise, rounds to up to 4 decimal places without trailing zeros (e.g. 8.8333)
    """
if abs(weight - round(weight)) < 1e-7: #check whether the weight is sufficiently close to an integer
  return str(int(round(weight))) #convert the rounded value to an integer and then to a string

# Round to 4 decimals and strip trailing zeros/dot
formatted = f"{weight:.4f}".rstrip('0').rstrip('.')
return formatted


def main():
  matrix = parse_input()
if not matrix:
  return

adj = upgma(matrix)

# Print adjacency list sorted by parent node ID and child node ID
for u in sorted(adj.keys()):
  adj[u].sort(key=lambda x: x[0]) #sort the children of each parent by their node ids to ensure that edges are printed in ascending childe-node order
for v, weight in adj[u]: #iterate through every child and its corresponding edge weight
  print(f"{u}->{v}:{format_weight(weight)}") #print required format: parent -> child : weight


if __name__ == "__main__":
  
#Exercise 2
  import sys


def parse_input():
  """Parses the distance matrix from standard input.
    
    Supports inputs separated by whitespace (tabs or spaces), with or
    without an explicit dimension n on the first line.
    """
raw_lines = [line.strip() for line in sys.stdin.read().strip().splitlines() if line.strip()]
if not raw_lines:
  return []

first_line_tokens = raw_lines[0].split() #split the first line into seperate token using whitespace

# If the first line contains a single integer, it is 'n'
if len(first_line_tokens) == 1:
  n = int(first_line_tokens[0])
matrix_lines = raw_lines[1:n + 1]
else:
  # Matrix starts directly on line 0
  n = len(raw_lines)
matrix_lines = raw_lines[:n]

matrix = []
for line in matrix_lines:
  row = [float(x) for x in line.split()]
matrix.append(row)

return matrix


def neighbor_joining_rec(D, clusters, next_node_id, adj):
  """Recursive helper for Neighbor Joining."""
n = len(D) #get number of active nodes remaining in current matrix

# Base case: 2 nodes remaining
if n == 2: #if only two nodes remain connect them directly
  c1, c2 = clusters[0], clusters[1] #obtain labels corresponding to the two remaining nodes
dist = D[0][1] #distance between the two remaining nodes is their connecting edge

# Directed edge from greater label to smaller label
parent = max(c1, c2)
child = min(c1, c2)

if parent not in adj:
  adj[parent] = []
adj[parent].append((child, dist))

return next_node_id

# 1. TotalDistance D(i)
total_dist = [sum(D[i]) for i in range(n)] #total_dist[i] represents the sum of distance from node i to every pther node in the current distance matrix

# 2. Construct Neighbor-Joining matrix D*
D_star = [[0.0] * n for _ in range(n)]
min_val = float("inf")
i_min, j_min = -1, -1

for i in range(n):
  for j in range(n):
  if i != j: #a node cannot be paired with itself
  D_star[i][j] = (n - 2) * D[i][j] - total_dist[i] - total_dist[j]
if D_star[i][j] < min_val: #check whether this value is smaller than the smallest value encountered so far
  min_val = D_star[i][j] #update smallest value
i_min, j_min = i, j #remember the indices of the pair producing this value

if i_min > j_min: #ensure that the smaller index comes first
  i_min, j_min = j_min, i_min

# 3. Compute Delta and Limb Lengths
delta = (total_dist[i_min] - total_dist[j_min]) / (n - 2)
limb_i = 0.5 * (D[i_min][j_min] + delta) #calculate limb length connecting node i to the new internal node
limb_j = 0.5 * (D[i_min][j_min] - delta) #calculate limb length connecting node j to the new internal node

# retrieve original labels of the two nodes being merged
node_i = clusters[i_min]
node_j = clusters[j_min]
node_m = next_node_id #assign an id to th enewly created internal node
next_node_id += 1

# 4. Construct new distance matrix
active_indices = [k for k in range(n) if k != i_min and k != j_min] #identify all nodes except the two selected nodes, which will remain active after the merge
m_distances = [] #list to store distances between each remaining node an newly created node m

for k in active_indices: #calculates the distance from every remaining noe k to the new node m
  d_km = 0.5 * (D[k][i_min] + D[k][j_min] - D[i_min][j_min])
m_distances.append(d_km)

new_D = []
for idx_r, r in enumerate(active_indices): #construct the rows corresponding to the remaining original nodes
  new_row = []
for c in active_indices: #copy the distances between the remaining original nodes
  new_row.append(D[r][c])
new_row.append(m_distances[idx_r]) #add distance from current node r to new node m
#idx_r identifies the corresponding distance in m_distances
new_D.append(new_row) #add completed row to the new distance matrix

m_row = m_distances + [0.0] #construct the final row corresponding to the new node m, its distance to the remaining nodes are stores in m_distances and its distance to itself is zero
new_D.append(m_row)

# 5. Update active clusters list by removing the two merged clusteres from active list and keeping all other clusters and append the new internal node m
new_clusters = [clusters[k] for k in active_indices] + [node_m]

# Recursive call
next_node_id = neighbor_joining_rec(new_D, new_clusters, next_node_id, adj)

# 6. Add directed limbs connecting greater label to smaller label
for child, weight in [(node_i, limb_i), (node_j, limb_j)]: #process both children and thei calculated limb lengths
  parent = max(node_m, child)
c = min(node_m, child)

if parent not in adj:
  adj[parent] = []
adj[parent].append((c, weight))

return next_node_id


def neighbor_joining(D):
  """Main Neighbor Joining function.

    Args:
        D: n x n distance matrix (list of lists).

    Returns:
        adj: Dict mapping parent_id -> list of (child_id, weight).
    """
n = len(D)
if n <= 1:
  return {}

clusters = list(range(n)) #initially every leaf represents its own cluster
next_node_id = n #internal node IDs begin at n, immediately after leaf ids
adj = {}

neighbor_joining_rec(D, clusters, next_node_id, adj)

return adj


# Alias for capital camel-case naming
NeighborJoining = neighbor_joining


def main():
  matrix = parse_input()
if not matrix:
  return

adj = neighbor_joining(matrix)

# Output directed adjacency list sorted by parent ID and child ID formatted to 3 decimal places
for u in sorted(adj.keys()):
  adj[u].sort(key=lambda x: x[0])
for v, weight in adj[u]:
  print(f"{u}->{v}:{weight:.3f}")


if __name__ == "__main__":
  main()