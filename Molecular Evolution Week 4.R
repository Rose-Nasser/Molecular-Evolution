#Week 4
#Exercise 1
import sys


# Standard integer mass table for the 20 amino acids
AMINO_ACID_MASSES = {
  'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
  'T': 101,'C': 103, 'I': 113, 'L': 113, 'N': 114,
  'D': 115,'K': 128, 'Q': 128, 'E': 129, 'M': 131,
  'H': 137,'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

def spectrum_graph(*args): #defines a function that constructs a spectrum graph
  # Handle various input formats (string, list, or multiple arguments)
  if len(args) == 1: #checks whether the function received exactly one argument
  arg = args[0] #stores the first and only argument in the variable arg
if isinstance(arg, str): #checks if the argument is a string
  nums = [int(x) for x in arg.strip().split()] #removes whitespeace from the beginning and end, splits the string into separate values, and conerts each value as an integer
elif isinstance(arg, (list, tuple)):
  nums = [int(x) for x in arg] #converts every element in the list or tuple into an integer
else: #if function received more than one argument
  nums = [int(arg)] #converts every argument into an integer and stores the results in a list
else:
  nums = [int(x) for x in args]

spectrum = [0] + sorted(nums) #adds 0 to the spectrum because the spectrum graph starts at mass 0
#sorted(nums) arranges all the input masses in ascending order

# Map each mass to a single representative amino acid to avoid duplicate edges
mass_to_aa = {}
for aa, mass in AMINO_ACID_MASSES.items():
  if mass not in mass_to_aa:
  mass_to_aa[mass] = aa  # Keeps the first encountered (e.g., 'K' for 128, 'I' for 113)

edges = [] #creates an empty list which will store all edges in the spectrum
n = len(spectrum) #calculates the total number of masses in spectrum
for i in range(n):
  for j in range(i + 1, n): #iterates through the masses that come after the current mass, prevents comparing a mass with itself and avoids generating duplicate edges in the opposite direction
  u = spectrum[i] #stores the current starting mass in variable v
v = spectrum[j] #stores the next mass in variable v
diff = v - u

if diff in mass_to_aa:
  edges.append((u, v, mass_to_aa[diff]))

return edges

#Exercise 2
import sys
from collections import deque

# Mass (int) to Amino Acid(s) (str)
AMINO_ACID_MASSES = {
  57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V',
  101: 'T', 103: 'C', 113: 'IL', 114: 'N', 115: 'D',
  128: 'KQ', 129: 'E', 131: 'M', 137: 'H', 
  147: 'F', 156: 'R', 163: 'Y', 186: 'W'
}

def compute_ideal_spectrum(peptide):
  # Ideal spectrum includes masses of all prefix and suffix subpeptides plus 0
  spectrum = {0} #creates a set containing 0, representing the starting mass, which automatically prevents duplicate mass values
n = len(peptide) #calcultes the total number of amino acids in the peptide

# Prefix masses
curr_mass = 0
for i in range(n):
  curr_mass += AA_TO_MASS.get(peptide[i], 0) #.e=get() returns 0 if the amino acid is not found in the dictionary
spectrum.add(curr_mass)

# Suffix masses
curr_mass = 0
for i in range(n - 1, -1, -1): #iterates backward through the peptide starting from last amino acid
  #n - 1 is the last index, and -1 means the loop stps before index -1
  curr_mass += AA_TO_MASS.get(peptide[i], 0)
spectrum.add(curr_mass)

return sorted(list(spectrum))

def decode_ideal_spectrum(*data):
  # .data() allows function to handle different data input formats (string, list, or multiple arguments)
  if len(data) == 1:
  item = data[0] #stores the only argument in the variable item
if isinstance(item, str):
  spectrum = [int(x) for x in item.strip().split()]
elif isinstance(item, (list, tuple)):
  spectrum = [int(x) for x in item]
else:
  spectrum = [int(item)]
else: #if function receives multiple separate arguments
  spectrum = [int(x) for x in data] #converts all arguments into integers and stores them in a list

target_spectrum = sorted(list(set([0] + spectrum))) #adds 0 to the input spectrum because every ideal spectrum starts at mass 0
# set() removes duplicate masses, list() converts the set back into a lits, sorted() arranges the masses in ascending order
max_mass = target_spectrum[-1] #gets the largest mass in the spectrum, this represnets the total mass of the peptide we want to reconstruct 

# Build adjacency list for the graph
graph = {node: [] for node in target_spectrum}
for i in range(len(target_spectrum)):
  for j in range(i + 1, len(target_spectrum)):
  diff = target_spectrum[j] - target_spectrum[i]
if diff in AMINO_ACID_MASSES:
  for aa in AMINO_ACID_MASSES[diff]:
  graph[target_spectrum[i]].append((target_spectrum[j], aa)) #adds a directed edge from the smaller mass to the larger mass. the edge is repesented as a tuple containing the destination mass and the amino acid labeling the edge

# Sort neighbors in descending order of target node mass to prioritize longer initial hops (e.g., 128 before 71)
for node in graph:
  graph[node].sort(key=lambda x: x[0], reverse=True) #sorts the outgoing edges according to their destiantion mass
#reverse = TRUE means the largest destination mass comes first

# BFS queue stores (current_node, current_peptide)
queue = deque([(0, "")]) #creates a BFS containing the initial state. The starting node is mass 0. The peptide is initially an empty string because no amino acids have been added

while queue: #contains searching as long as the queue contains unexplored states
  current_node, current_peptide = queue.popleft()
#current_node is the current spectrum mass
#current_peptide is the amino acid sequence constructed so far
#popleft() implements the FIFO behavior of BFS
if current_node == max_mass: #checks whether the current path has reached the largest spectrum mass
  if compute_ideal_spectrum(current_peptide) == target_spectrum: #calculates the ideal spectrum of the reconstructed peptide, compares it with the original target spectrum. the peptide is accepted only if both spectra are identical 
  return current_peptide
continue #if the maximum mass was reached but the spectrum does not match, this path is rejected and the algorithm moves to the next queued state

for next_node, aa in graph.get(current_node, []): #gets all outgoing edges from the current node, otherwise returns an empty list
  queue.append((next_node, current_peptide + aa)) #adds a new state to the end of the queue. the new state contains the destination mass and the peptide formed by adding the amino acid to the existing sequence

return None

#Exercise 3
import sys


# Standard integer mass table for the 20 amino acids (plus toy alphabet X=4, Z=5)
AMINO_ACID_MASSES = {
  'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
  'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115,
  'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137,
  'F': 147, 'R': 156, 'Y': 163, 'W': 186,
  'X': 4, 'Z': 5
}

# Insert your peptide_to_vector function here, along with any subroutines you need
def peptide_to_vector(peptide): #converts peptide sequence into a binary peptide vector
  #input is a string containing amino acid letters and output is a lits of 0s and 1s representing the peptide's prefix masses
  
  # Compute prefix masses
  prefix_masses = set() #empty set to store cumulative masses of peptide prefixes, and automatically removes duplicate values
current_mass = 0
for amino_acid in peptide:
  current_mass += AMINO_ACID_MASSES[amino_acid]
prefix_masses.add(current_mass)

total_mass = current_mass

# Construct the binary vector from mass 1 to total_mass
peptide_vector = [1 if m in prefix_masses else 0 for m in range(1, total_mass + 1)]
 #creates a binary vector containing one position for every integer mass from 1 up to and including the total peptide mass
return peptide_vector

#Exercise 4
import sys
# Please do not remove package declarations because these are used by the autograder.

# Mass-to-amino-acid lookup (for ambiguous masses 113 and 128, any choice is accepted)
MASS_TO_AMINO_ACID = {
  4: 'X', 5: 'Z', 57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V',
  101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D',
  128: 'Q', 129: 'E', 131: 'M', 137: 'H',
  147: 'F', 156: 'R', 163: 'Y', 186: 'W'
}

def vector_to_peptide(peptide_vector):
  # Identify all prefix masses (1-indexed positions where value is 1)
  prefix_masses = [i + 1 for i, val in enumerate(peptide_vector) if val == 1] #identifies all prfix masses (1- indexed positions where value is 1)
#enumarete peptide_vector returns both the index and value of each element
#i represents the index starting from 0 and val represents the value at that index
# the resulting list contains all prefix masses

# Reconstruct amino acids from consecutive prefix mass differences
peptide = []
prev_mass = 0
for m in prefix_masses:
  amino_acid_mass = m - prev_mass #calculates the mass difference between the current suffix mass and the previous prefix mass. This difference represents the mass of one amino acid
peptide.append(MASS_TO_AMINO_ACID.get(amino_acid_mass, ''))
prev_mass = m

return "".join(peptide)

#Exercise 5
import sys
# Please do not remove package declarations because these are used by the autograder.

# Standard integer mass table for the 20 amino acids
STANDARD_AMINO_ACID_MASSES = {
  'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
  'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115,
  'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137,
  'F': 147, 'R': 156, 'Y': 163, 'W': 186
}


def run_dp(spectrum, m, masses): #finds highest scoring peptide using a dynamic programming
  #spectrum: the spectral vector, including the initial zero
  #m: the total mass, equal to number of spectral vectors positions
  #masses: dictionary containing amino acids and their integer masses
  dist = [-float('inf')] * (m + 1) #creates a list of length m+1
#dist[v] stores the max score achievable when reaching mass v
#initially, all masses are assigned negative infinity, meaning that they have not yet been reached
parent = [-1] * (m + 1) #creates a list that stores the previous mass used to reach each mass
#this allows us to reconstruct the optimal peptide later
#-1 means that no previous mass has been recorded
backtrack_aa = [''] * (m + 1) #creates a list that stores the amino acid used to reach each mass, which also helps reconstruct the the peptide after completing the DP algorithm


dist[0] = 0

for u in range(m + 1): #iterates through all possible mass positions, from 0 to m, u represents the current mass position
  if dist[u] == -float('inf'): #checks whether the current mass is unreachable
  continue
for aa, mass in masses.items():
  v = u + mass
if v <= m: #checks that new mass does not exceed max mass
  if dist[u] + spectrum[v] > dist[v]: #calculates the score of the new path
  #dist[u] is the best score up to the current mass
  #spectrum[v] is the score contributed by reaching mass v
  #checks whether this path produces a better score than the score already stored for mass v
  dist[v] = dist[u] + spectrum[v] #updates max score for reaching mass v
parent[v] = u #records the previous mass u
backtrack_aa[v] = aa #records the aa that connects u to v

if dist[m] == -float('inf'): #checks whether the final mass m is unreachable
  return None

curr = m
peptide = []
while curr > 0:
  peptide.append(backtrack_aa[curr])
curr = parent[curr]

peptide.reverse()
return "".join(peptide)

def peptide_sequencing(*args): #*args allows different input fromats
  # Handle different possible argument formats from the autograder safely
  if len(args) == 1:
  arg = args[0]
if isinstance(arg, str):
  spectrum_str = arg
elif isinstance(arg, (list, tuple)):
  spectrum_str = " ".join(map(str, arg)) #converts every eleemnt into a string and joins the elemetns using spaces
else: #if neither string, list nor tuple
  spectrum_str = str(arg)
else: #if receives multiple arguments
  spectrum_str = " ".join(map(str, args))

# Parse the spectral vector and prepend s0 = 0
spectrum = [0] + list(map(int, spectrum_str.strip().split()))
#removes surrounding whitespace from input string and splits tsring into individual values. Converts each alue into an integer
#adds 0 at the beginning because spectral vectors start at mass position 0

m = len(spectrum) - 1 #calculates the max mass position. Since initial zero was added, we substract one from the list length

# Try with standard amino acid masses first
result = run_dp(spectrum, m, STANDARD_AMINO_ACID_MASSES)
if result is not None:
  return result

# Fall back to extended masses if standard ones cannot form a path
return run_dp(spectrum, m, EXTENDED_AMINO_ACID_MASSES)

#Exercise 6
import sys

# Please do not remove package declarations because these are used by the autograder.

# Insert your suffix_array function here, along with any subroutines you need
def suffix_array(text):
  """
    Constructs the suffix array for a given string Text.
    Returns a list of starting indices of the sorted suffixes.
    """
# Create a list of tuples containing each suffix and its original starting index
suffixes = [(text[i:], i) for i in range(len(text))]

# Sort suffixes lexicographically
suffixes.sort(key=lambda x: x[0])

# Extract and return the indices
return [index for suffix, index in suffixes]

