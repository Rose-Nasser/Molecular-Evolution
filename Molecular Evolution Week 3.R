#Week 3 
#Exercise 1
import sys
import itertools
import re
from collections import defaultdict


def small_parsimony(leaf_strings, adjacency_list):
  alphabet = "ACGT"
n = len(leaf_strings) #counts number of leaf nodes (the nodes whose DNA sequences are already known)

children = defaultdict(list)
nodes = set(range(n)) #creates a set containing all leaf node numbers: 0, 1, ..., n-1
child_nodes = set()

# Parse the adjacency list.
for edge in adjacency_list:
  if isinstance(edge, str): #checks if edge is provided as a string, such as "4 -> 0"
  u, v = edge.split("->", 1) #splits the string into 2 parts: parent node u and child node v
else:
  u, v = edge #assigns the first element to the parent and second element to the child

u, v = int(u), int(v)
children[u].append(v) #adds node v to list of children belonging to parent u
child_nodes.add(v) #records node v as a child node
nodes.update((u, v)) #adds both nodes to set of all nodes in the tree

roots = nodes - child_nodes #finds root by subtracting all child nodes from complete set of nodes. The root is never a child
if len(roots) != 1: #check if there is exactly one root in the tree
  raise ValueError("Expected exactly one root.")

root = next(iter(roots)) #extracts the only root node from set

if not leaf_strings: #checks whether the input contains any leaf DNA sequences
  raise ValueError("No leaf DNA strings provided.")

if len(set(map(len, leaf_strings))) != 1: #calculates length of every DNA sequence and checks whether all lengths are identical
  raise ValueError("All leaf DNA strings must have equal lengths.")

# Postorder traversal.
order = []

def postorder(v): #a function that visits a node's children before visiting the node itself
  for child in children[v]:
  postorder(child)
order.append(v) #adds current node after all its children have been visited

postorder(root)

if set(order) != nodes: #check if every node in the tree was visited during traversal
  raise ValueError("The tree contains disconnected nodes.")

for v in nodes:
  if v < n: #checks whether the node is a leaf (leaf nodes are numbered from 0 to n-1)
  if children[v]: #checks if leaf noe has any children
  raise ValueError("Leaf nodes must be numbered 0 to n-1.")
elif len(children[v]) != 2: #check if an internal node has exactly 2 children
  raise ValueError("Each internal node must have two children.")

labels = {v: [] for v in nodes} # creates a dictionary where each node will store its assigned nucleotide at every DNA position
total_score = 0
m = len(leaf_strings[0]) #gets the number of nucleotide positions in each DNA sequence

# Solve each DNA position independently.
for pos in range(m):
  scores = {} #dictionary to score minimum parsimony score for each node and each possible nucleotide
backtrack = {} #dictionary to remember which nucleotides were selected for the children of each internal node

for v in order:
  scores[v] = {}

if v < n: #check if current node is a leaf
  observed = leaf_strings[v][pos] #retrieves the known nucleotide at the current DNA position for this leaf
for k in alphabet:
  scores[v][k] = 0 if k == observed else float("inf") #assigns a score of 0 to observed nucleotide and infinity to all other nucelotides
else: #if current node is internal node
  left, right = children[v] #retrieves the two children of current internal node
backtrack[v] = {} #dictionary to store best nucleotide choices for both children

for k in alphabet: #calculate cost for every possinle nucelotide i
  left_cost = { #dictionary containing possible costs for assigning each nucleotide to the left child
    i: scores[left][i] + (i != k)
    for i in alphabet
  }
right_cost = { #creates a dictionary containing the possible costs for assigning each nucleotide to the right child
  j: scores[right][j] + (j != k) #adds 1 if the right child's nucleotide differs from the parent's nucleotide , otherwise adds 0
  for j in alphabet
}

scores[v][k] = ( #stores minimum total score for current node assuming its nucleotide is k
  min(left_cost.values()) #find minimum possible value for left subtree
  + min(right_cost.values()) #adds minimum possible score for right subtree
)

backtrack[v][k] = ( #records the best nucelotide assignments for the 2 children
  min(alphabet, key=left_cost.get), #finds the nucelotide that gives the lowest score for the left child
  min(alphabet, key=right_cost.get) #finds the nucelotide that gives the lowest score for the right child
)

root_char = min(alphabet, key=lambda k: scores[root][k]) #selects the nucelotide that gives the minimum parsimony score at the root
total_score += scores[root][root_char] #adds minimum score for this DNA position to the overall parsimony score

# Backtrack to assign characters to every node.
assigned = {root: root_char} #assigns the selected nucleotide to the root node
stack = [root] #creates a stack containing the root, which will be used to traverse the tree

while stack: #continues untill all nodes have received their nucelotide assignemnts
  v = stack.pop() #removes the most recently added node from the stack

if v in backtrack: #checks whether the current node is an internal node with recorded child assignemnts 
  left, right = children[v] #retrieves the left and right children of the current node
assigned[left], assigned[right] = ( #assigns the optimal nucelotides to both children
  backtrack[v][assigned[v]] #retrieves the two best child nucelotides corresponding to the current parent's assigned nucleotide 
)
stack.extend([left, right]) #adds both children to the stack so that their descendants can be processed

for v in nodes:
  labels[v].append(assigned[v]) #adds the nucelotide assigned to node v at the current DNA position

# Return the score and labels separately, as the tester expects.
labels = {
  v: "".join(chars) #joins all nucleotides assigned to node v into one DNA string
  for v, chars in labels.items() #repeats this operation for every node in the labels dictionary 
}

ordered_labels = [labels[v] for v in sorted(labels)] #creates a list of DNA sequences ordered by their node numbers

return int(total_score), ordered_labels


#Exercise 2
import sys
from collections import defaultdict
import queue
import itertools
import re


def small_parsimony_unrooted_tree(leaf_strings, adjacency_list):
  alphabet = "ACGT"
n = len(leaf_strings)

if n == 0: #check whether the input contains any DNA sequences
  raise ValueError("No DNA strings provided.")

if len(set(map(len, leaf_strings))) != 1: #calculates the length of every DNA sequence and checks whether all sequences have the same length
  raise ValueError("All DNA strings must have equal lengths.")

# Construct an undirected tree.
graph = defaultdict(set) #creates an adjacency dictioanry where each node stores a set of its neighboring nodes. Sets prevent duplicate neighbors
nodes = set(range(n))

for edge in adjacency_list:
  if isinstance(edge, str):
  edge = edge.strip()
if not edge or "->" not in edge: #checks if edge is empty or if it has no "->" separator
  continue
u, v = edge.split("->", 1) #splits the edge into 2 parts: u is the first node and v is the second node
u, v = int(u.strip()), int(v.strip()) #removes surrounding spaces and converts both node identifiers into integers
else: #if edge is already represented as a pair such as (4, 0)
  u, v = int(edge[0]), int(edge[1]) #extracts both node identifiers and converts them into integers

graph[u].add(v)
graph[v].add(u)
nodes.update((u, v))

# An unrooted binary tree has internal nodes numbered n or higher.
internal_nodes = sorted(v for v in nodes if v >= n) #identifies internal nodes. leaves are numbered 0 to n-1, so internal nodes have IDs >= n. they are sorted in increasing order

if not internal_nodes:
  raise ValueError("No internal node found.")

# Temporarily root the tree at an internal node.
root = internal_nodes[0] #selects the internal node with the smallest ID as a temporary root
children = defaultdict(list) #creates a dictionary that will store the children of each node after temporarily rooting the tree
parent = {root: None} #records the parent of each node. The temporary node has no parent
preorder = []
stack = [root]

while stack:
  v = stack.pop()
preorder.append(v)

for w in sorted(graph[v], reverse=True): #iterates through the neighbors of v in descending order to make traversal deterministic 
  if w == parent[v]: #check whether the neighbor is the parent of the current node
  continue

if w in parent: #checks whether the neighbor has already been assigned a parent
  raise ValueError("The adjacency list contains a cycle.")

parent[w] = v #recirds v as the parent of neighbor w
children[v].append(w)
stack.append(w)

if set(preorder) != nodes: #compare visited nodes will all nodes in the tree
  raise ValueError("The tree is disconnected.")

# Reverse preorder to obtain postorder.
order = list(reversed(preorder)) #reverses preorder ro obtain a child-before-parent processing order for this rooted tree

# Validate leaves and internal nodes.
for v in nodes:
  if v < n:
  if len(graph[v]) != 1: #checks if leaf has exactly one neighbor in an unrooted tree
  raise ValueError(f"Leaf {v} must have degree one.")
elif len(graph[v]) != 3: #checks whether an internal node has exactly 3 neighbors
  raise ValueError(
    f"Internal node {v} must have degree three."
  )

sequence_length = len(leaf_strings[0]) #gets the number of nucleotide positions in each DNA sequence
labels = {v: [] for v in nodes} #creates a dictionary where each node will store its assigned nucleotide at every DNA position
total_score = 0

# Apply Small Parsimony independently at each DNA position.
for pos in range(sequence_length):
  scores = {}
backtrack = {}

for v in order:
  scores[v] = {}

if v < n:
  # Leaf: its nucleotide is fixed.
  observed = leaf_strings[v][pos] # retrieves the known nucleotide of leaf v at the current DNA position

for k in alphabet:
  scores[v][k] = 0 if k == observed else float("inf") #gives the observed nucleotide a score of zero and all other nucleotides an infinite score because the leaf's nucleotide cannot change 

else:
  backtrack[v] = {}

# The temporary root has three children;
# all other internal nodes have two.
for k in alphabet:
  scores[v][k] = 0 #initializes teh total score for assigning nucleotide k to node v
backtrack[v][k] = [] #there are 2 children for ordinary internal nodes and 3 for the temporary root

for child in children[v]:
  costs = {
    i: scores[child][i] + (i != k) #adds the child's subtree score and a penalty of 1 if the child's nucleotide differs from parent's nucleotide
    for i in alphabet
  }

best_char = min(alphabet, key=lambda i: costs[i]) #find the nucleotide that gives the lowest total cost for this child
scores[v][k] += costs[best_char] #adds the minimum code of this child to the total score of the current parent node
backtrack[v][k].append(best_char) #saves the best nucleotide  assigned to this child so it can be reconstructed later

# Select an optimal nucleotide for the temporary root.
root_char = min(alphabet, key=lambda k: scores[root][k]) #selects the nucleotide that gives the minimum parsimony score at the temporary root
total_score += scores[root][root_char] #adds the minimum score for the current DNA position to the overall parsimony score

# Backtrack to reconstruct every node's nucleotide.
assigned = {root: root_char} #assigns the selected nucleotide to the temporary root
stack = [root]

while stack:
  v = stack.pop()

if v not in backtrack: #checks if current node is a leaf. Leaves do not have recorded child assignemnts 
  continue

for child, char in zip( #pairs each child with its previously calculated optimal nucleotide
  children[v], backtrack[v][assigned[v]] #retrieves the children of v and the best nucleotide choices corresponding to the assigned nucleotide of v
):
  assigned[child] = char #assigns the selected nucleotide to the child node
stack.append(child)

for v in nodes:
  labels[v].append(assigned[v]) #adds the nucleotide assigned to node v at the current DNA position

# Return labels in increasing node-ID order.
ordered_labels = [
  "".join(labels[v]) for v in sorted(nodes)
]

return int(total_score), ordered_labels


# Optional compatibility alias.
def small_parsimony(leaf_strings, adjacency_list):
  return small_parsimony_unrooted_tree(
    leaf_strings, adjacency_list
  )

#Exercise 3
import sys
from collections import defaultdict
import queue
import itertools
import re


def nearest_neighbors(a, b, adjacency_list):
  a, b = int(a), int(b)

# Build an ordered, undirected adjacency dictionary that will represent the tree. Each node will be associated with a list of its neighboring nodes
graph = {}

def ensure_node(u): #ensures node exists in the graph
  if u not in graph:
  graph[u] = []

def add_edge(u, v): #adds an undirected edge between nodes u and 
  u, v = int(u), int(v)
ensure_node(u)
ensure_node(v) #ensure both nodes exist in the graph

if v not in graph[u]: #checks whether v is already listed as a neighbor of u
  graph[u].append(v) #adds v to the neighbor list of u if it is not already present
if u not in graph[v]: #checks whether u is already listed as a neighbor of v
  graph[v].append(u) #adds u to the neighbor list of , making the edge undirected

def parse_adjacency(obj): #defines a recursive helper function that reads different possible formats of adjacency-list input
  flat_numbers = [] #to store node numbers when they are provided as consecutive pairs

if isinstance(obj, dict): #checks whether adjacency list is provided as a dictionary
  for u, neighbors in obj.items(): #iterates through each dictionary entry, where u is a node and neighbors contains its adjacent nodes
  u = int(u) #converts the current node identifier into an integer
ensure_node(u) #makes sure the current node exists in the graph

if isinstance(neighbors, str): #checks whether neighbors are provided as a string
  matches = re.findall( #searches string for directed edge representations
    r'(\d+)\s*->\s*(\d+)', neighbors #matches 2 groups of digits separated by an arrow, allowing optional spaces
  )
if matches:#checks whether any complete edges were found in the string
  for x, y in matches: 
  add_edge(int(x), int(y)) #adds each extracted edge to the undirected graph
else: #if no complete arrow-formatted edge was found
  for v in re.findall(r'\d+', neighbors): #extracts every sequence of digits from the neighbor string 
  add_edge(u, int(v)) #adds each extracted number as a neighbor of the current node u

elif isinstance(neighbors, int): #checks whether the neighbor is provided as a single integer
  add_edge(u, neighbors) #adds that integer as a neighbor of u

else: #if neighbors is another iterable, sich as a list or tuple
  for v in neighbors:
  if isinstance(v, str) and "->" in v:
  match = re.fullmatch( #checks whether the entire string follows the expected edge format
    r'\s*(\d+)\s*->\s*(\d+)\s*', v
  )
if match:
  add_edge(
    int(match.group(1)),
    int(match.group(2)) #extracts and converts the first and second node identifier
  )
else: #executed if the neighbor is not an arrow-formatted string
  add_edge(u, int(v)) #adds v as a neighbor of the current node u

return

if isinstance(obj, str): #checks whether the adjacency list is provided as a single string
  matches = re.findall(r'(\d+)\s*->\s*(\d+)', obj) #extracts all edges represented using the arrow notation

if matches:
  for u, v in matches:
  add_edge(int(u), int(v)) #adds each edge to the undirected graph
else: #if no arrow-formatted edges were found
  flat_numbers.extend( #adds all extracted numbers to the temporary list
    int(x) for x in re.findall(r'\d+', obj) #extracts every sequence of digits and converts it into an integer
  )

elif isinstance(obj, (list, tuple)):
  for item in obj:
  if isinstance(item, str):
  matches = re.findall( #searches the string for an arrow-formatted edges
    r'(\d+)\s*->\s*(\d+)', item #identifies two node numbers separated by an arrow
  )
if matches: #checks whether any edges were found
  for u, v in matches:
  add_edge(int(u), int(v))
elif item.strip().isdigit(): #checks whether the string contains only a node number, ignoring surrounding spaces
  flat_numbers.append(int(item.strip())) #converts the number into an integer and stores it temporarily

elif isinstance(item, int):
  flat_numbers.append(item)

elif isinstance(item, (list, tuple)) and len(item) == 2: #checks whether the element is a pair containing exactly two node identifiers 
  add_edge(int(item[0]), int(item[1])) #adds the pair as an edge between its two nodes

elif isinstance(item, dict): #checks whether element is another dictionary
  parse_adjacency(item)

if flat_numbers: #checks whether any node numbers were collected without explicit edge separators
  if len(flat_numbers) % 2 != 0: #checks whether the number of collected integers is odd, meaning one node has no matching partner
  raise ValueError("Invalid adjacency list: unmatched node ID.")

for i in range(0, len(flat_numbers), 2): #iterates through the list  two numbers at a time
  add_edge(flat_numbers[i], flat_numbers[i + 1]) #treats 2 consecutive pair as an edge and adds it to the graph

parse_adjacency(adjacency_list)

# Check the specified central edge.
if a not in graph or b not in graph[a]: #checks whether node a exists and whether b is one of its neighbors
  raise ValueError(
    f"The specified edge ({a}, {b}) was not found in the tree."
  )

if len(graph[a]) != 3 or len(graph[b]) != 3: #checks whether both endpoints of the selected edge have exactly 3 neighbors
  raise ValueError(
    "Both endpoints of the specified edge must have degree 3."
  )

# Find the two branches attached to each endpoint, excluding a-b.
a_branches = [v for v in graph[a] if v != b] #retrieves the 2 neighbors of node a other than b. These represent 2 branches attached to a
b_branches = [v for v in graph[b] if v != a] #retrieves the 2 neighbors of node b a other than a. These represent 2 branches attached to b

# Keep a_branches[0] fixed; swap a_branches[1] with each
# of the two branches adjacent to b.
branch_a = a_branches[1] 
results = []

for branch_b in b_branches:#selects one of the 2 branches attached to node a. This branch will be used in both rearrangements
  new_graph = { #creates a copy of the original graph so that the original tree remains unchanged
    u: list(neighbors) #copies the neighbor list of each node into a new independent list
    for u, neighbors in graph.items()b #repeats this operation for every node in the original graph
  }

# Replace the branch at each endpoint and its reciprocal edge.
replacements = [ #creates a list describing the 4 neighbor-list changes required for the branch interchange
  (a, branch_a, branch_b), #at node a, replaces branch_a with branch_b
  (branch_a, a, b), #at branch_a, replaces its connection to a with a connection to b
  (b, branch_b, branch_a), #at node b, replaces branch_b with branch_a
  (branch_b, b, a) #at branch_b, replaces its connection to b with a connection to a
]

for node, old_neighbor, new_neighbor in replacements: #iterates through all 4 changes needed to preserve the undirected edges
  idx = new_graph[node].index(old_neighbor) #finds the position of the old neighbor in the current node's neighbor list
new_graph[node][idx] = new_neighbor #replaces the old neighbor with the new neighbor at the same position

results.append(new_graph) #adds the newly generated tree to the results list


return results[0], results[1] #returns the 2 alternative adjacency dictionaries generated by the nearest-neighbor interchange


# Compatibility alias, if needed.
def nearest_neighbors_of_a_tree(a, b, adjacency_list):
  return nearest_neighbors(a, b, adjacency_list)

#Exercise 4
import sys
import re
import itertools
import queue
from collections import defaultdict
from collections.abc import Mapping


def _normalize_node(x): #converts node identifiers into a consistent format
  if isinstance(x, int):
  return x
x = str(x).strip() #conerts x into a string and removes unnecessary spaces around it 
if re.fullmatch(r"\d+", x): #checks whether the entire string consists of one or more digits
  return int(x)
return x


def _parse_adjacency(adj): #conerts different adjacency list formats into an ordered undirected graph
  """Convert supported adjacency formats into an ordered undirected graph."""
graph = {} #dictionary representing the tree, each noe will be associated with a list of its neighbors

def ensure(u): #ensure node u exists in the graph
  if u not in graph: #if node u has already been added
  graph[u] = [] #creates an empty neighbor list for node u if it does not exist 

def add_edge(u, v): #adds an undirected edge between two nodes
  u, v = _normalize_node(u), _normalize_node(v) #normalize both node identifiers so numeric strings become integers and DNA strings become strings
if u == v: #checks whether an edge connects to a node itself
  raise ValueError("Self-loops are not allowed.")
ensure(u)
ensure(v)
if v not in graph[u]: #check if v is already a neighbor of u
  graph[u].append(v) #adds v to the neighbor list of u if it is not already present
if u not in graph[v]: #checks whether u is already a neighbor of v
  graph[v].append(u) #adds u to the neighbor list of v, making the edge undirected

def parse(obj): #defines a recursive function that processes different possible formats of adjacency-list input
  if isinstance(obj, Mapping): #checks whether the input is dictionary like
  for u, neighbors in obj.items(): #iterates through each dictionary entry, where u is  anode and neighbors represents its connected nodes
  u = _normalize_node(u) #normalize the current node identifier
ensure(u) #makes sure the current node exists in the graph

if isinstance(neighbors, str): #checks whether neighbor information is provided as a string
  matches = re.findall( #searches the string for edges written using arrow notation
    r"([A-Za-z]+|\d+)\s*->\s*([A-Za-z]+|\d+)", #matches alphabetic or numeric node identifiers separeted by an arrow, allowing optional spaces
    neighbors #specifies the string in which the pattern is searched
  )
if matches: #if arrow formatted edges were found
  for x, y in matches: #iterates through every extracted edge
  add_edge(x, y) #adds the extracted edge to the undirected graph
else: #executes this block if no arroww formatted edge was found
  add_edge(u, neighbors) #treats the entire neighbor string as a single neighbor of u

elif isinstance(neighbors, int):
  add_edge(u, neighbors)

else: #if the neighbors are provided as an iterable, such as a list or tuple
  for v in neighbors:
  if isinstance(v, str) and "->" in v:
  parse(v)
else:
  add_edge(u, v)
return

if isinstance(obj, str):
  matches = re.findall( #searches the string for all arrow formatted edges
    r"([A-Za-z]+|\d+)\s*->\s*([A-Za-z]+|\d+)", #defines the pattern for identifying two node identifiers separated by an arrow
    obj #specifies the input string to search
  )
for u, v in matches:
  add_edge(u, v)
return

if isinstance(obj, (list, tuple)):
  if (len(obj) == 2 #checks whetehr list or tuple contains exactly two elements
      and all(isinstance(x, (str, int)) for x in obj) #checks whether both elements are strings or integers
      and not any("->" in str(x) for x in obj)): #ensure that neither element contains an arrow, meaning the pair can represent one edge
  add_edge(obj[0], obj[1]) #treats the two elements as the endpoints of an edge and adds that edge
else: #if input is not a simple two-element edge pair
  for item in obj: #iterates through every element in the list or tuple
  parse(item) #recursively processes each element according to its data type
return

raise ValueError("Unsupported adjacency-list format.")

parse(adj) #parses the adjacency list supplied to the function

if not graph: #checks if parsing process produced an empty graph
  raise ValueError("The adjacency list is empty.")

return graph


def _small_parsimony_score(graph): #defines a function thta calculates the minimum parsimony score and reconstructs the DNA sequence of every node
  """Return minimum score and a mapping from node IDs to DNA labels."""
alphabet = "ACGT"
nodes = list(graph) #creates a list containing all nodes in graph
leaves = [v for v in nodes if isinstance(v, str)] #identifies leaf nodes
internal = sorted(v for v in nodes if isinstance(v, int)) #identifies internal nodes. there identifiers are integers sorted in ascending order

if not leaves or not internal: #checks if tree contains both leaves and internal nodes
  raise ValueError("The tree must contain leaves and internal nodes.")

if len({len(v) for v in leaves}) != 1: #calculates the ;engths of all leaf DNA sequences and checks whether they  are identical
  raise ValueError("All leaf DNA strings must have equal length.")

for v in nodes: #iterates through every node in the tree
  expected = 1 if v in leaves else 3 #assigns the expected degree: leaves must have one neighbor, while internal nodes must have three
if len(graph[v]) != expected: #compares the actual number of neighbors with the expected degree
  raise ValueError(
    f"Node {v} has degree {len(graph[v])}; expected {expected}." #provides the node identifier, its actual degree and its expected degree
  )

root = internal[0] #selects the internal node with the smallest identifier as a temporary root

# Root the unrooted tree temporarily.
parent = {root: None} #creates a dictionary that records the parent of each node. the temporary root has no parent
children = {v: [] for v in nodes} #dictionary where every node initially has an empty list of children
preorder = [] #creates an empty list to store nodes in preorder traversal
stack = [root] #initailizes a stack containing the temporary root

while stack:
  v = stack.pop()
preorder.append(v)

for w in graph[v]: #iterates through all neighbors of current node
  if w == parent[v]: #checks whether the neighbor is the parent of the current node
  continue
if w in parent: #checks whether the neighbor has already been visited and assigned a parent
  raise ValueError("The adjacency list contains a cycle.")

parent[w] = v
children[v].append(w)
stack.append(w)

if len(preorder) != len(nodes): #compares the number of visited nodes with the total number of nodes
  raise ValueError("The tree is disconnected.")

postorder = list(reversed(preorder)) #reverses preorder to obtain an order in which descendants are processed before their ancestors
labels = {v: [] for v in nodes} #creates a dictionary where each node will store its assigned nucleotide at every DNA position
total_score = 0
sequence_length = len(leaves[0]) #gets the number of nucelotide positions in each leaf DNA sequence

# Dynamic programming independently at every DNA position.
for pos in range(sequence_length): #repeats the scoring calculation independently for every nucelotide position
  scores = {} #dictionary to score minimum parsimony score for each node and each possible nucleotide

for v in postorder: #processes nodes in postorder, ensuring children are evaluated before their parent
  scores[v] = {} #creates an empty dictionary to store the 4 nucelotide scores for node v

if isinstance(v, str): #checks whether current node is a leaf
  observed = v[pos] #retrieves the observed nucelotide at the current position from the leaf's DNA sequence
for k in alphabet:
  scores[v][k] = ( #assigns the score for each possible nucelotide at this leaf
    0 if observed == k else float("inf") #gives the observed nuceltoide a score of 0 and all other nucelotides an infinite score
  )
else: #if current node is an internal node
  for k in alphabet: #calculates the minimum score assuming the current internal node has nucelotide k
  scores[v][k] = sum( #calculates the total score by adding the minimum contributions from all children
    min( #finds minimum possible cost for one child
      scores[ch][i] + (i != k) #adds the child's subtree score and a penalty of 1 if the child's nucelotide differs from parent's nucelotide
      for i in alphabet
    )
    for ch in children[v]
  )

root_char = min(alphabet, key=lambda k: scores[root][k]) #finds the nucelotide that gives the minimum parsimony score at the temporary root
total_score += scores[root][root_char] #adds the minimum score for the current DNA position to the overall parsimony score

# Trace back the optimal assignment to every node.
assigned = {root: root_char} #assigns the selected nucleotide to the temporary root
stack = [root] #initializes a stack with the root to begin the backtracking process

while stack: #continue untill all nodes have received their nucleotide assignment 
  v = stack.pop() #removes the most recently added node from stack
 
for ch in children[v]: #iterate through all children of the current node
  assigned[ch] = min( #assigns the nucelotide that minimizes the child's score given the parent's assigned nucleotide
    alphabet, #consider all possible nucleotides
    key=lambda i: scores[ch][i] + (i != assigned[v]) #selects the nucelotide with the lowest subtree score + the penalty for differing from the parent
  )
stack.append(ch)

for v in nodes:
  labels[v].append(assigned[v]) #adds the nucelotide assigned to node v at the current DNA position

labels = {v: "".join(chars) for v, chars in labels.items()} #joins the nucelotide characters collected at all positions to reconstruct each node's complete DNA sequence
return int(total_score), labels #returns the total minimum parsimony score and a dictionary mapping every node to its reconstructed DNA sequence


def _nearest_neighbor_trees(graph, a, b): #generates the two alternative trees obtained by applying NNI to the internal edge connecting a and b
  """Generate the two nearest-neighbor trees around internal edge a-b."""
if b not in graph[a]: #checks whether b is actually connected to a
  return [] #returns an empty list if the specified edge does not exist

a_branches = [v for v in graph[a] if v != b] #identifies the two branches connected to a, excluding the central edge a-b
b_branches = [v for v in graph[b] if v != a] #identifies the two branches connected to b, excluding the central edge a-b

if len(a_branches) != 2 or len(b_branches) != 2:#checks whether both endpoints have exactly 2 additional branches
  return []

moving_branch = a_branches[1] #selects the second branch attached to a. This branch will be exchanged with each branch attached to b
neighbors = [] #empty list to store the 2 aletrnative trees

for branch_b in b_branches: #iterates through the 2 branches attached to b, producing one alternative tree fro each
  new_graph = { #creates a copy of the original graph so that the original tree is not modified
    u: list(neighbors) #copies the neighbor list of each node into a new independent list
    for u, neighbors in graph.items() #repeats this operation for every node in the original graph
  }

# Swap edge a-moving_branch with edge b-branch_b.
new_graph[a][new_graph[a].index(moving_branch)] = branch_b #at node a, replaces the connection to moving_branch with a connection to branch_b
new_graph[moving_branch][new_graph[moving_branch].index(a)] = b #at moving_branch, replaces its connection to a with a connection to b

new_graph[b][new_graph[b].index(branch_b)] = moving_branch #at node b, replaces the connection to branch_b with a connection to moving_branch
new_graph[branch_b][new_graph[branch_b].index(b)] = a #at branch_b, replaces its connection to b with a connection to a

neighbors.append(new_graph) #adds newly generated tree to list o flaternative trees

return neighbors


def large_parsimony_nni(n, adj): #searches for a more parsimonious tree using nearesr-neighbor interchange
  """
    Return a list of (score, labels, adjacency) triples.
    Includes the initial tree and every strictly improving NNI step.
    """
graph = _parse_adjacency(adj) #converts the supplied adjacency list into an undirected graph

leaves = [v for v in graph if isinstance(v, str)] #identifies the leaf nodes, which are represented by DNA strings
if len(leaves) != int(n): #checks whetehr the number of identifies leaves matches the expected number n
  raise ValueError(
    f"Expected {n} leaves but found {len(leaves)}."
  )

score, labels = _small_parsimony_score(graph) #calculates the initial's tree minimum parsimony score and reconstructs the DNA labels of all nodes

# The testing harness requires triples: (score, labels, adjacency).
steps = [] #empty list to store the improving trees and their corresponding results

while True: #repeats the search til no neighboring tree has a strictly lower parsimony score
  best_score = score
best_graph = None

internal_nodes = sorted(
  v for v in graph if isinstance(v, int) #selects an internal node identifiers because internal nodes are represented by integers
)

# Examine every internal edge and its two NNI neighbors.
for a in internal_nodes: #iterates through every internal node as the first endpoint of a possible central edge
  for b in graph[a]: #iterates through all neighbors of node a
  if not isinstance(b, int) or a >= b: #skips leaf neighbors and prevents examining the same undirected edge twice
  continue 

for candidate in _nearest_neighbor_trees(graph, a, b): #generates the two NNI alternative trees arounf the internal edge a-b
  candidate_score, _ = _small_parsimony_score(candidate) #calculates the parsimony score of the candidate tree
#the underscore means the reconstructed labels are not needed at this stage

if candidate_score < best_score:
  best_score = candidate_score
best_graph = candidate

# Stop when no neighbor strictly improves the score.
if best_graph is None:
  break

graph = best_graph
score, labels = _small_parsimony_score(graph) #recalculates the score and reconstructs all node labels for the newly selected tree

steps.append( #adds the current improving tree and its information to the results list
  (
    score, #stores the minimum parsimony score of the current tree
    dict(labels), #stores a copy of the dictionary containing the reconstructed DNA sequences for all nodes
    {u: list(vs) for u, vs in graph.items()} #stores a copy of the adjacency dictionary representing the current tree
  )
)

return steps
