# Project Report - Alignment

## Baseline

### Design Experience

I did this with Spencer Zaugg on 10/27/2025.

My plan is to use an np.array to store the memoized data of the edit distance for different substrings. I will initialize it by filling in the base cases and side base cases. For back pointers, I will use a hash map to point each entry ij of the matrix to its previous entry that it descended from. I will keep this map updated as I fill in entries.

After finding the bottom right corner of the matrix. I will then recurse through the back pointer dictionary to assemble the reverse order of operations. After that, I will iterate over those in reverse order and apply the necessary adjustments to both strings depending on if it was diagonal, down, right. This should make it O(mn).

He wanted to use a numpy array.

### Theoretical Analysis - Unrestricted Alignment

#### Time 

##### align - **O(mn)**

For my analysis, m will be the length of the first string, and n will be the length of the second.

```py
calculator = AlignmentCalculator(seq1, seq2, banded_width, # O(m)
                                indel_penalty, match_award, 
                                sub_penalty, gap_open_penalty)
dist, prev = calculator.calculate_matrix() # O(mn) see below
edits = get_edit_sequence(seq1, seq2, prev) # O(m + n) see below
aligned_seq1, aligned_seq2 = get_sequence_alignment(seq1, seq2, gap, edits) # O(m + n)
score = math.inf if aligned_seq1 is None else dist[len(seq1)-1][len(seq2)-1] # O(1)
return score, aligned_seq1, aligned_seq2
```

I have not include the parameters and the rest of the method signature because its really big. This align function instantiates an AlignmentCalculator object (I did this to avoid code duplication and have better design) and uses it to calculate both the distance adjacency list and the back pointer prev adjacency list.

Both of these operations are O(m) and O(mn) respectively (see my analysis of those functions below). Getting the sequence of edits and reconstructing the alignment sequence are both O(m+n) (see below).

This means that the AlignmentCalculator.calculate_matrix() function defines our big O at **O(mn)**

```py
class AlignmentCalculator:
    def __init__(self, seq1, seq2, 
                 banded_width, indel_penalty, 
                 match_award, sub_penalty, gap_open_penalty):
        self.seq1 = seq1
        self.seq2 = seq2
        self.banded_width = banded_width
        self.indel_penalty = indel_penalty
        self.match_award = match_award
        self.sub_penalty = sub_penalty
        self.gap_open_penalty = gap_open_penalty
        self.dist = {i: {} for i in range(-1, len(self.seq1))} # O(m)
        self.prev = {}
        self.mode = {i : {} for i in range(-1, len(self.seq1))} # O(m)
```

For initialization, we simply assign references to input values. However, it does m operations to initialize our distance and mode matrices.

```py
def calculate_matrix(self):
        for i, j in self.position_iterator(): # O(mn)
            self.calculate_distance(i, j) # O(1)
        return self.dist, self.prev

def position_iterator(self):
    if self.banded_width == -1:
        for i in range(-1, len(self.seq1)):    #O(mn)
            for j in range(-1, len(self.seq2)):
                yield i, j
    else:
        for i in range(-1, len(self.seq1)):    # Explained in Core
            j_start = max(-1, i - self.banded_width)
            j_end = min(len(self.seq2) - 1, i + self.banded_width)
            for j in range(j_start, j_end + 1):
                yield i, j~
```

Here, we see that our matrix calculation function simply loops over all m * n positions in the dynamic programming matrix. This means that this function will scale by at least O(mn). See my justification of the calculate_distance function as being constant time below.

```py
def calculate_distance(self, i, j):
    if i == -1 and j == -1: # O(1)
        self.dist[i][j] = 0
        self.mode[i][j] = "m"
    elif i == -1:               # O(1)
        gap_modifier = self.gap_open_penalty if j == 0 else 0
        self.dist[i][j] = self.dist[i][j-1] + self.indel_penalty + gap_modifier
        self.prev[i, j] = (i, j-1)
        self.mode[i][j] = "i"
    elif j == -1:     # O(1)
        gap_modifier = self.gap_open_penalty if i == 0 else 0
        self.dist[i][j] = self.dist[i-1][j] + self.indel_penalty + gap_modifier
        self.prev[i, j] = (i-1, j)
        self.mode[i][j] = "d"
    else:
        self.get_best_distance(i, j) # O(1)
```

In this functions, the operations that are taking place are all constant time. We have several comparisons and floating point arithmetic, so we consider those to be constant. Reading and accessing a hashmap data structure is also a constant time operation. I have justified the get_best_distance method below.

```py
def get_best_distance(self, i, j):
    left_scenario, diag_scenario, min_scenario = self.get_scenarios(i, j) # O(1)
    if min_scenario == diag_scenario: # O(1)
        self.prev[i, j] = (i-1, j-1)
        self.mode[i][j] = "m"
    elif min_scenario == left_scenario: # O(1)
        self.prev[i, j] = (i, j-1)
        self.mode[i][j] = "i"
    else: # O(1)
        self.prev[i, j] = (i-1, j)
        self.mode[i][j] = "d"
```

Similar to the other example. This function only performs hash table lookups and floating point arithmetic operations. This functions remains constant time. See get_scenarios below.

```py
def get_scenarios(self, i, j):
    if j in self.dist[i-1]: # O(1)
        gap_modifier = self.gap_open_penalty if self.mode[i-1][j] != "d" else 0
        up_scenario = self.dist[i-1][j] + self.indel_penalty + gap_modifier
    else: # O(1)
        up_scenario = math.inf
    if j-1 in self.dist[i]: # O(1)
        gap_modifier = self.gap_open_penalty if self.mode[i][j-1] != "i" else 0
        left_scenario = self.dist[i][j-1] + self.indel_penalty + gap_modifier
    else: # O(1)
        left_scenario = math.inf
    diag_modifier = self.match_award if self.seq1[i] == self.seq2[j] else self.sub_penalty
    diag_scenario = self.dist[i-1][j-1] + diag_modifier # O(1)
    min_scenario = min(up_scenario, left_scenario, diag_scenario) # O(1)
    self.dist[i][j] = min_scenario # O(1)
    return left_scenario, diag_scenario, min_scenario
```

This only thing noteworthy here is calculating the minimum, which will always be of only three values. The rest of the operations are simple hash table lookups and variable reference reassignments. This make the function constant time. Because all the constituent function that calculate_distance depends on run in constant time (and it only performs constant time operations using these functions), it means that calculate_distance must be constant time as well.

```py
def get_edit_sequence(seq1, seq2, prev):
    rev_edits = []
    i, j = len(seq1) - 1, len(seq2) - 1 # O(1)
    while True:  # O(m + n) see below
        if not prev.get((i, j)): # O(1)
            break
        i_prev, j_prev = prev[(i, j)]
        if i_prev == i - 1 and j_prev == j - 1:  # O(1)
            rev_edits.append("sub")
        elif i_prev == i and j_prev == j - 1:  # O(1)
            rev_edits.append("ins")
        else:                                 # O(1)
            rev_edits.append("del")
        i, j = i_prev, j_prev
    edits = rev_edits[::-1] #O(m + n)
    return edits
```

We now examine the time complexity for computing the sequence of edits prior to assemble the sequence itself. As we recurse (logically, not by functional recursion) through the prev back pointer hash map, each recursion through is a constant time operation. This is because its simply hash map lookups and edits.

However, we have to understand that the maximum amount of times we iterate through the while loop is bounded by the maximum number of edits. This must then logically be the height of our matrix m + the width of the matrix n. This is because m + n is the worst case scenario for alignment. This makes this function bound by O(m + n)

```py
def get_sequence_alignment(seq1, seq2, gap, edits):
    if (not edits) and (len(seq1) != 0 and len(seq2) != 0):
        return None, None
    aligned_seq1, aligned_seq2 = "", ""
    seq1_shifts, seq2_shifts = 0, 0
    for edit in edits: # O(m + n)
        if edit == "sub":
            aligned_seq1 += seq1[seq1_shifts]
            aligned_seq2 += seq2[seq2_shifts]
            seq1_shifts += 1
            seq2_shifts += 1
        elif edit == "ins":
            aligned_seq1 += gap
            aligned_seq2 += seq2[seq2_shifts]
            seq2_shifts += 1
        else:
            aligned_seq1 += seq1[seq1_shifts]
            aligned_seq2 += gap
            seq1_shifts += 1
    return aligned_seq1,aligned_seq2
```

Within each iteration, the operations are trivial constant time hash table lookups and increments. This means that this function is bound by the number of edits that have happened. We established previously that the maximum edits it m + n. That makes this function O(m + n)

Between all the high level functions and operations of this algorithm, its clear that the calculate_matrix() function defines the big O for this algorithm. That means that algin has a time complexity of **O(mn)**

#### Space

##### align - **O(mn)**

```py
calculator = AlignmentCalculator(seq1, seq2, banded_width, # see below
                                indel_penalty, match_award, 
                                sub_penalty, gap_open_penalty)
dist, prev = calculator.calculate_matrix() # O(mn)
edits = get_edit_sequence(seq1, seq2, prev) # O(mn)
aligned_seq1, aligned_seq2 = get_sequence_alignment(seq1, seq2, gap, edits) # O(m + n)
score = math.inf if aligned_seq1 is None else dist[len(seq1)-1][len(seq2)-1]
return score, aligned_seq1, aligned_seq2
```

The distance and previous matrices take up O(mn) space because they represent m x n matrices. The aligned sequences themselves are bound by the total number of edits manuevers in the algorithm, which is O(m + n). I will examine the AlignmentCalculator space separately.

```py
class AlignmentCalculator:
    def __init__(self, seq1, seq2, 
                 banded_width, indel_penalty, 
                 match_award, sub_penalty, gap_open_penalty):
        self.seq1 = seq1
        self.seq2 = seq2
        self.banded_width = banded_width
        self.indel_penalty = indel_penalty
        self.match_award = match_award
        self.sub_penalty = sub_penalty
        self.gap_open_penalty = gap_open_penalty
        self.dist = {i: {} for i in range(-1, len(self.seq1))} # O(mn)
        self.prev = {}                                         # O(mn)
        self.mode = {i : {} for i in range(-1, len(self.seq1))}
```

We see that the AlignmentCalculator has the dist, prev (I will discuss mode in core, but just know that ir grows to mn too), which take up O(mn) space. This is what eventually gets returned.

```py
def calculate_matrix(self):
        for i, j in self.position_iterator(): #O(1)
            self.calculate_distance(i, j)
        return self.dist, self.prev

def position_iterator(self):
    if self.banded_width == -1:
        for i in range(-1, len(self.seq1)):
            for j in range(-1, len(self.seq2)):
                yield i, j                        #O(1)
    else:
        for i in range(-1, len(self.seq1)):
            j_start = max(-1, i - self.banded_width)
            j_end = min(len(self.seq2) - 1, i + self.banded_width)
            for j in range(j_start, j_end + 1):
                yield i, j                        #O(1)
```

As a construction itself, the iterator takes constant space because it doesn't store all the values in any data structure: it simply evaluates the next iteration when asked. In a for loop, the previous values are garbage collected by the python virtual machine. This makes these functions constant in space.

```py
def calculate_distance(self, i, j): # O(1) see below
    if i == -1 and j == -1:
        self.dist[i][j] = 0
        self.mode[i][j] = "m"
    elif i == -1:
        gap_modifier = self.gap_open_penalty if j == 0 else 0
        self.dist[i][j] = self.dist[i][j-1] + self.indel_penalty + gap_modifier
        self.prev[i, j] = (i, j-1)
        self.mode[i][j] = "i"
    elif j == -1:
        gap_modifier = self.gap_open_penalty if i == 0 else 0
        self.dist[i][j] = self.dist[i-1][j] + self.indel_penalty + gap_modifier
        self.prev[i, j] = (i-1, j)
        self.mode[i][j] = "d"
    else:
        self.get_best_distance(i, j)
```

This only modifies member data structures, which we have already established are O(mn) in space. All the intermediary values for computation here are fixed and do not scale. They are therefore constant in space.

```py
def get_best_distance(self, i, j): # O(1)
    left_scenario, diag_scenario, min_scenario = self.get_scenarios(i, j)
    if min_scenario == diag_scenario:
        self.prev[i, j] = (i-1, j-1)
        self.mode[i][j] = "m"
    elif min_scenario == left_scenario:
        self.prev[i, j] = (i, j-1)
        self.mode[i][j] = "i"
    else:
        self.prev[i, j] = (i-1, j)
        self.mode[i][j] = "d"
```

This is the exact same. It simply store a few intermediary values but nothing that scales with our variables, so this is constant.

```py
def get_scenarios(self, i, j): # O(1)
    if j in self.dist[i-1]:
        gap_modifier = self.gap_open_penalty if self.mode[i-1][j] != "d" else 0
        up_scenario = self.dist[i-1][j] + self.indel_penalty + gap_modifier
    else:
        up_scenario = math.inf
    if j-1 in self.dist[i]:
        gap_modifier = self.gap_open_penalty if self.mode[i][j-1] != "i" else 0
        left_scenario = self.dist[i][j-1] + self.indel_penalty + gap_modifier
    else:
        left_scenario = math.inf
    diag_modifier = self.match_award if self.seq1[i] == self.seq2[j] else self.sub_penalty
    diag_scenario = self.dist[i-1][j-1] + diag_modifier
    min_scenario = min(up_scenario, left_scenario, diag_scenario)
    self.dist[i][j] = min_scenario
    return left_scenario, diag_scenario, min_scenario
```

As you can see, this is merely an extension of what we have seen before. The member data structure are modified up to O(mn) space, but this function itself only needs a fixed number of intermediary values for the computation. This is also O(1) space.

```py
def get_sequence_alignment(seq1, seq2, gap, edits): # O(m + n)
    if (not edits) and (len(seq1) != 0 and len(seq2) != 0):
        return None, None
    aligned_seq1, aligned_seq2 = "", ""
    seq1_shifts, seq2_shifts = 0, 0
    for edit in edits:  # (m + n)
        if edit == "sub":
            aligned_seq1 += seq1[seq1_shifts]
            aligned_seq2 += seq2[seq2_shifts]
            seq1_shifts += 1
            seq2_shifts += 1
        elif edit == "ins":
            aligned_seq1 += gap
            aligned_seq2 += seq2[seq2_shifts]
            seq2_shifts += 1
        else:
            aligned_seq1 += seq1[seq1_shifts]
            aligned_seq2 += gap
            seq1_shifts += 1
    return aligned_seq1,aligned_seq2  # O(2m + 2n)
```

We have a fixed set of values used for intermediary computation that are garbage collected at each iteration of the for loop. However, we are storing both the edits, of which there are O(m + n) maximum, and the aligned sequences, which would then be O(2m + 2n) at the maximum.

Because of this intermediary storing of the edits and sequences, this function has a space complexity of O(m + n)

```py
def get_edit_sequence(seq1, seq2, prev):
    rev_edits = []
    i, j = len(seq1) - 1, len(seq2) - 1
    while True:
        if not prev.get((i, j)):
            break
        i_prev, j_prev = prev[(i, j)]
        if i_prev == i - 1 and j_prev == j - 1:
            rev_edits.append("sub")
        elif i_prev == i and j_prev == j - 1:
            rev_edits.append("ins")
        else:
            rev_edits.append("del")
        i, j = i_prev, j_prev
    edits = rev_edits[::-1] # O(m + n)
    return edits
```

This is largely repetition of what I said before. The reversal of the edits and the edits themselves both take up O(m + n) space because that is the maximum number of edits.

To conclude, we see that the adjacency list data structures grow the most in the AlignmentCalculator. This then binds the Big O space complexity for the whole align function. I conclude that this would then be **O(mn)** space.

### Empirical Data - Unrestricted Alignment

| N    | time (ms) |
|------|-----------|
| 500  |           |
| 1000 |           |
| 1500 |           |
| 2000 |           |
| 2500 |           |
| 3000 |           |


### Comparison of Theoretical and Empirical Results - Unrestricted Alignment

- Theoretical order of growth: 
- Empirical order of growth (if different from theoretical): 


![](fill-me-in.png)

*Fill me in*

## Core

### Design Experience

I did this with Spencer Zaugg on 10/27/2025.

Instead of using an mxn matrix to store the information. I will use a hashmap instead. I anticipate that the best thing to do would be to have the keys and values be integers to correspond to the index of each respective string in the theoretical matrix.

However, in my iteration, I will not iterate over all points. I will make my iteration smart enough to only iterate over the banded section and only consider the two(in the edge cases) side available to it.

I will still maintain a hashmap of back pointers as usual, but I will need to reformat the dictionary keys and values with some delimiter (likely a tuple) to make sure that no data is lost. I will then use this hashmap to construct the alignments.

He also wanted to use a dictionary.

### Theoretical Analysis - Banded Alignment

#### Time 

*Fill me in*

#### Space

*Fill me in*

### Empirical Data - Banded Alignment

| N     | time (ms) |
|-------|-----------|
| 100   |           |
| 1000  |           |
| 5000  |           |
| 10000 |           |
| 15000 |           |
| 20000 |           |
| 25000 |           |
| 30000 |           |

### Comparison of Theoretical and Empirical Results - Banded Alignment

- Theoretical order of growth: 
- Empirical order of growth (if different from theoretical): 


![](fill-me-in.png)

*Fill me in*

### Relative Performance Of Unrestricted Alignment versus Banded Alignment

*Fill me in*


## Stretch 1

### Design Experience

I did this with Spencer Zaugg on 10/27/2025.

I will simply parse the fasta file by using string manipulation techniques available in python such as split (and regular expressions if necessary). After that, I will simply run the edit distance algorithm, comparing unknown to each of the sequences to find which one is the closest. I will print out the sequences and their alignment.

### Code

```python
# Fill me in
```

### Alignment Scores

*Fill me in*

## Stretch 2

### Design Experience

I did this with Spencer Zaugg on 10/27/2025.

This will be quite similar to other implementations, but I will use another matrix (or hashmap or np.array) to determine whether or the candidate previous box in question was inside of a first gap, second gap, match (non-penalized substitution), or penalized substitution. If so, it will inform the cost (via the gap_open_penalty) of taking that route.

This will inform me regarding the state of the analysis. I will then adjust the prospective score for these various points on that basis. This will occur for each box.

### Alignment Outcome Comparisons

##### Sequences and Alignments

*Fill me in*

##### Chosen Parameters and Better Alignments Discussion

*Fill me in*

## Project Review

*Fill me in*
