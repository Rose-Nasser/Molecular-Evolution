#Week 5
#Exercise 1
import sys
# Please do not remove package declarations because these are used by the autograder.

# Standard integer mass table for the 20 amino acids (plus toy X=4, Z=5)
AMINO_ACID_MASSES = {
  'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
  'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115,
  'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137,
  'F': 147, 'R': 156, 'Y': 163, 'W': 186,
  'X': 4, 'Z': 5
}

def peptide_identification(spectrum, proteome): #defines a fucntion that identifies the highest scoring peptide substring within given proteome, using a spectral vector
  best_score = float('-inf') #initialize the highest score to negative infinity. This enusres that the score of any valid peptide will be higher initially
best_peptide = "" #initialize an empty string to store the best peptide found so far

target_mass = len(spectrum) #determines the required total mass of the peptide. the length of the spectral vector represents the target mass

for start in range(len(proteome)): #Iterates through every possible starting position in proteome
  #resets the accumulated mass to zero for the new starting position for a peptidesubstring
  current_mass = 0 # resets the accumulated mass to 0 fro the new starting position
current_score = 0 #resets the accumulated score to 0 for the new pwptide substring

for end in range(start, len(proteome)): #extends the peptide one amino acid at a time, starting from current start position
  aa = proteome[end] #retrieves the amino acid at the current end position
current_mass += AMINO_ACID_MASSES[aa] #looks up the mass of the current amino acid in the dictionary and adds it to the total mass of the current peptide

if current_mass > target_mass:
  break #if the accumulated mass exceeds the required target mass

current_score += spectrum[current_mass - 1] #adds the spectral score corresponding to the current cumulative mass
#substracting 1 converts the mass into a zero-based python list index

# Only consider peptides with exactly the required mass
if current_mass == target_mass:
  if current_score > best_score:
  best_score = current_score
best_peptide = proteome[start:end + 1] #extracts and stores the peptide substring from the proteome
#+1 is used to include the amino acid at position end

break

return best_peptide

#Exercise 2
import sys

def psm_search(SpectralVectors, Proteome, threshold):
  
  mass = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97,
    'V': 99, 'T': 101, 'C': 103, 'I': 113,
    'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137,
    'F': 147, 'R': 156, 'Y': 163, 'W': 186,
    'X': 4, 'Z': 5
  }


if isinstance(SpectralVectors, str):
  SpectralVectors = list(map(int, SpectralVectors.split()))
elif SpectralVectors and isinstance(SpectralVectors[0], str): #checks if Spectralvectors is a non-empty list whose first element is a string
  SpectralVectors = [
    list(map(int, v.split())) for v in SpectralVectors
  ]
#processes each spectral vector separately
#v represents one spectral vector in string format
#v.split() separates its values
#map(int, ...) converts the values into integers

if SpectralVectors and isinstance(SpectralVectors[0], int): #checks if SpectralVectors is a non empty list containing integers rather than a list of lists. This means that only one spectral vector was provided 
  SpectralVectors = [SpectralVectors] #wraps the single spectral vector inside another list
#This gives the input a consistent structure: multiple vectors are represented as a list of lists. even one vector is now treated as a list containing one vector

result = set()

for vector in SpectralVectors:
  
  target_mass = len(vector) 
best_score = float('-inf')
best_peptide = None

for start in range(len(Proteome)):
  
  current_mass = 0
score = 0

for end in range(start, len(Proteome)):
  
  aa = Proteome[end] #retrieves the aa located at the current end position
current_mass += mass[aa] #looks up the mass of the current aa in the mass dictionary and adds this mass to the total mass of the current peptide

if current_mass > target_mass:
  break

score += vector[current_mass - 1] #adds the spectral score corresponding to the current cumulative mass

if current_mass == target_mass:
  
  if score > best_score:
  best_score = score
best_peptide = Proteome[start:end + 1] #extracts and stores the current peptide substring

break

if best_peptide is not None and best_score >= threshold:
  result.add(best_peptide)

return result


PSMSearch = psm_search

#Exercise 3
import sys

# Standard integer mass table (plus toy X=4, Z=5)
AMINO_ACID_MASSES = {
  'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
  'T': 101, 'C': 103, 'I': 113, 'L': 113,
  'N': 114, 'D': 115, 'K': 128, 'Q': 128,
  'E': 129, 'M': 131, 'H': 137, 'F': 147,
  'R': 156, 'Y': 163, 'W': 186,
  'X': 4, 'Z': 5
}


def spectral_alignment(peptide: str, spectral_vector: list[int], k: int) -> str:
  """
    Find the highest-scoring peptide related to the original
    peptide by at most k modifications.

    A modification changes an amino acid's mass by a nonzero
    integer delta. Every modified mass must be positive.
 """
#k is the max number of allowedaa modifications
#-> str means that the function returns a string


n = len(spectral_vector)
L = len(peptide)

if L == 0 or L > n:
  return ""

original_masses = []

for aa in peptide:
  if aa not in AMINO_ACID_MASSES:
  raise ValueError(f"Unknown amino acid: {aa}")
original_masses.append(AMINO_ACID_MASSES[aa]) #retrieves the original mass of the amino acid and adds it to the list

k = min(k, L) #limits the maximum number of modifications to the peptide length. each aa can be modified at most once

NEG_INF = float('-inf')

# dp[j][t] = maximum score after processing the current
# number of amino acids, with total mass j and t modifications.
dp = [[NEG_INF] * (k + 1) for _ in range(n + 1)] #creates a 2D dynamic programming table
#rows represent the accumulated peptide mass (0 to n)
#columns represent the number of modifications (0 to k)
dp[0][0] = 0 #total mass, number of modifications and accumulated score are all 0

# Store predecessors for reconstructing the answer.
parent = {}

for i in range(L):
  original_mass = original_masses[i] #retrieves the original mass of aa at position i

new_dp = [[NEG_INF] * (k + 1) for _ in range(n + 1)] #creates a new DP table for the next iteration
#stores the best scores after processing one additional aa
#using a separate table prevents the same aa from being processed multiple times within the same iteration

for j in range(n + 1):
  for t in range(k + 1):
  
  current_score = dp[j][t] #retrieves the best score for the current state. the state represents accumulated mass j and t modifications

if current_score == NEG_INF:
  continue

# At least one unit of mass is needed
# for every remaining amino acid.
remaining = L - i - 1 #calculates how many amino acids remain after the current one 
max_mass = n - j - remaining #calculates the  max mass that can be assigned to the current aa
#It reverses at least  one mass unit for every remaining aa so that the final total mass does not exceed n

if max_mass < 1:
  continue

# Try every possible positive mass.
for new_mass in range(1, max_mass + 1):
  
  modified = (new_mass != original_mass) #checks whether the new mass differes from original mass

new_t = t + int(modified) #updates the number of modifications. The modifications count increase by one only when the aa mass changes

if new_t > k:
  continue

new_j = j + new_mass #calculates the new accumulated peptide mass after adding the current amino acid

# Spectral vector uses 1-based prefix masses.
new_score = current_score + spectral_vector[new_j - 1]
#adds the spectral score sorresponding to the new cumulative mass
#new_j - 1 converts the cumulative mass into a zero-based python list index

if new_score > new_dp[new_j][new_t]: #checks whether this new path produces a higher score
  
  new_dp[new_j][new_t] = new_score #updates the DP table witht the higher score

parent[(i + 1, new_j, new_t)] = (
  j, t, new_mass
)
#stores the previous state and the mass selected for the current aa
dp = new_dp

# Find the best score with total mass exactly n
# and at most k modifications.
best_score = NEG_INF
best_t = -1

for t in range(k + 1):
  if dp[n][t] > best_score: #looks at the DP score for total mass n and t modifications
  #checks whether it is better than the best score found so far
  best_score = dp[n][t] #updates the highest final spectral score
best_t = t

if best_t == -1:
  return ""

# Backtrack to reconstruct the modified peptide.
result = []
j = n
t = best_t

for i in range(L, 0, -1): #moves backwards through the peptide
  
  previous_j, previous_t, new_mass = parent[(i, j, t)]
#retrieves the previous DP state and the mass assigned to the aa at the current position
original_mass = original_masses[i - 1] #retrieves the original mass of the aa
delta = new_mass - original_mass

aa = peptide[i - 1] #retrieves the original aa letter at the current position

if delta != 0:
  result.append(f"{aa}({delta:+d})") #stores the aa together with its mass modification
else:
  result.append(aa)

j = previous_j
t = previous_t

result.reverse()

return ''.join(result)

#Exercise 4
import sys


def spectral_dict_prob(spectral_vector: list[int], threshold: int,
                       max_score: int) -> float:
  #defines a function that calculates the probability of a spectral dictionary
  # -> float indicates that the function returns a decimal probability
  
  n = len(spectral_vector)

if max_score < 0:
  return 0.0

# Select the amino acid alphabet.
if n < 57:
  amino_acid_masses = [4, 5] #uses the toy amino acid alphabet, X has a mass of 4 and Z has a mass of 5
aa_probability = 0.5 
else: #if target mass is at least 57 use the standard 20 amino acid alphabet
  amino_acid_masses = [
    57, 71, 87, 97, 99, 101, 103,
    113, 113, 114, 115, 128, 128,
    129, 131, 137, 147, 156, 163, 186
  ]
aa_probability = 1.0 / 20

# dp[mass][score] = probability
# Only scores from 0 to max_score are stored.
dp = [
  [0.0] * (max_score + 1)
  for _ in range(n + 1)
]
#creates a 2D dynamic programming table
#rows repr4esent the accumulated peptide mass (0 to n)
#columns represent the accumulated spectral score (0 to max_score)
#every cell is initially set to 0.0 because no probabilities have been calculated yet


# Empty peptide.
dp[0][0] = 1.0 #total mass and total score euqal to 0 and probability = 1 (100%)

# Dynamic programming.
for mass in range(n + 1):
  
  for score in range(max_score + 1):
  
  probability = dp[mass][score] #retrieves the probability of reaching current state, which is defined by accumulated mass and accumulated score

if probability == 0.0:
  continue

for aa_mass in amino_acid_masses:
  
  new_mass = mass + aa_mass

if new_mass > n: #checks if the new peptide mass exceeds target mass. if it does, the peptide is discarded because it cannot return to a smaller mass
  continue

# Add the spectral score at the new prefix mass.
new_score = score + spectral_vector[new_mass - 1] #retrieves the spectral score corresponding to the new cumulative mass and adds it to the previously accumulated score

# Discard scores outside the allowed table height.
if 0 <= new_score <= max_score:
  
  dp[new_mass][new_score] += (
    probability * aa_probability
  )
#updates the probability of reachinf the new state
#probability = probability of reaching the previous state
#aa_probability = probability of selecting the new amino acid
#multiplying these values gives the probability of reaching the new state through this particular amino acid
# += adds this probability to the existing probability because multiple diffeent peptide sequences may reach the same mass and score

# Sum probabilities of all peptides meeting the threshold.
if threshold > max_score:
  return 0.0

answer = sum(dp[n][max(0, threshold):])
#calculates the final probability by adding the probabilities of all states that have the required total mass m
#dp[n] represents all possible scores for pepties with total mass n
#max(0, threshold) determines the first score index to include
#: means that all scores from this index to max_score are included
#therefore only peptides with total mass = n AND spectral score >== threshold contribute to the final answer
return answer

#exercise 5
import sys

# Standard integer mass table (plus toy X=4, Z=5)
AMINO_ACID_MASSES = {
  'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
  'T': 101, 'C': 103, 'I': 113, 'L': 113,
  'N': 114, 'D': 115, 'K': 128, 'Q': 128,
  'E': 129, 'M': 131, 'H': 137, 'F': 147,
  'R': 156, 'Y': 163, 'W': 186,
  'X': 4, 'Z': 5
}


def spectral_alignment(peptide: str, spectral_vector: list[int], k: int) -> str:
  """
    Find the highest-scoring peptide related to the original
    peptide by at most k modifications.

    A modification changes an amino acid's mass by a nonzero
    integer delta. Every modified mass must be positive.
    """
#peptide is the original amino acid sequence
#spectral_vector is a list containing the spectral scores
#k is the max number of allowed amino acid modifications
# -> str indicates that the function return a string

n = len(spectral_vector)
L = len(peptide)

if L == 0 or L > n:
  return ""
#checks whether the epptide is empty or has more amino acids than th etarget mass can accommodate 
#since every amino acid must have a posiitve mass, either condition makes a valid alignment impossible
original_masses = []

for aa in peptide:
  if aa not in AMINO_ACID_MASSES:
  raise ValueError(f"Unknown amino acid: {aa}")
original_masses.append(AMINO_ACID_MASSES[aa])

k = min(k, L) #makes usre that the max number of modifications does not exceed the number of amino acids
#each amino acid can be modified at most once. Therefore, k cannot meaningfully exceed L


NEG_INF = float('-inf')

# dp[j][t] = maximum score after processing the current
# number of amino acids, with total mass j and t modifications.
dp = [[NEG_INF] * (k + 1) for _ in range(n + 1)]
#creates a 2D dynamic programming table
#rows represent the accumulated peptide mass (0 to n)
#columns represent the number of modifications (0 to k
#each cell stores the max spectral score possible for that particular mass and modification count)
dp[0][0] = 0
#accumulated mass, number of modifications and accumulated spectral score all = 0

# Store predecessors for reconstructing the answer.
parent = {}

for i in range(L):
  original_mass = original_masses[i]

new_dp = [[NEG_INF] * (k + 1) for _ in range(n + 1)]
#creates a new DP table for the next iteration
#it will store the best scores after processing the current amino acid
# a separate table prevents the same amino acid from being processed multiple times in one iteration

for j in range(n + 1):
  for t in range(k + 1):
  
  current_score = dp[j][t] #retrieves the max score stored in the current DP state
#this state has accumulated mass j and t modifications

if current_score == NEG_INF:
  continue
#checks whether the current state is reachable
#if its score is negative inf, no valid path has reached this state


# At least one unit of mass is needed
# for every remaining amino acid.
remaining = L - i - 1 #calculates the number of amino acids remaining after processing the current aa
max_mass = n - j - remaining #calculates the max possible mass that can be assigned to the current aa
#n - j represents the mass still available
#remaining reverses at least one mass unit for every aa that has not yet been processed
#this prevents the algorithm from selecting a mass that makes it impossible to complete the peptide

if max_mass < 1:#checks whether there is enough mass available to assign a positive mass to the current amino acid
  continue

# Try every possible positive mass.
for new_mass in range(1, max_mass + 1): #tests every possible positive mass for the current amino acid
  #the new mass can be = or different from its original mass
  
  modified = (new_mass != original_mass) #returns TRUE if amino acid is modified and returns FALSE if the original mass is retained

new_t = t + int(modified)

if new_t > k:
  continue

new_j = j + new_mass

# Spectral vector uses 1-based prefix masses.
new_score = current_score + spectral_vector[new_j - 1]
#calculates the new accumulated spectral score
#current_score is the score accumulated previously 
#spectral_vector[new_j - 1] is the score corresponding to the new cumulative mass
if new_score > new_dp[new_j][new_t]: #checks whether the new path produces a higher score than the score currently stored in this DP cell
  
  new_dp[new_j][new_t] = new_score #updates the DP table with the higher score than the score than the score currently stored in this DP call
#udates the DP table with the higher score
#this ensures that each state stores its max possible spectral score
parent[(i + 1, new_j, new_t)] = (
  j, t, new_mass
)
#stores info about previous DP states
#i + 1 is the number of aa processed
#new_j is the new accumulated mass
#new_t is the new number of modifications
#j is a previous accumulated mass
#t is the previous number of modifications
#new_mass is the mass assigned to the current aa

dp = new_dp
#replaces the old DP table with the newly calculated table
#this prepares the algorithm to process the next aa
#after the final interaction, dp contains the best scores after processing all L aa

# Find the best score with total mass exactly n
# and at most k modifications.
best_score = NEG_INF #initializes the highest final spectral score to negative infinity
best_t = -1 

for t in range(k + 1):
  if dp[n][t] > best_score: #checks if the score for the target mass n and modification count t is higher than the best score so far
  best_score = dp[n][t] 
best_t = t

if best_t == -1: #checks if a valid final state was found
  return ""

# Backtrack to reconstruct the modified peptide.
result = []
j = n
t = best_t

for i in range(L, 0, -1): #iterates backwards through the aa
  #starts at L, meaning all aa have been processed
  #moves backward until reaching the first aa
  
  previous_j, previous_t, new_mass = parent[(i, j, t)]
#retrieves the info stored in the parent dictionary
#previous_j = accumulated mass before the current aa
#previous_t = number of modifications before the current aa
#this allows the algorithm to identify how the optimal peptide was constructed

original_mass = original_masses[i - 1] #retrieves the original mass of the current aa
delta = new_mass - original_mass

aa = peptide[i - 1]

if delta != 0: #checks if mass of aa was modified
  result.append(f"{aa}({delta:+d})")
#adds the aa and its mass modification to the result list
else:
  result.append(aa)

j = previous_j
t = previous_t

result.reverse()

return ''.join(result)
