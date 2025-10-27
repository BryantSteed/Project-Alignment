# Project Report - Alignment

## Baseline

### Design Experience

I did this with Spencer Zaugg on 10/27/2025.

My plan is to use an np.array to store the memoized data of the edit distance for different substrings. I will initialize it by filling in the base cases and side base cases. For back pointers, I will use a hash map to point each entry ij of the matrix to its previous entry that it descended from. I will keep this map updated as I fill in entries.

After finding the bottom right corner of the matrix. I will then recurse through the back pointer dictionary to assemble the reverse order of operations. After that, I will iterate over those in reverse order and apply the necessary adjustments to both strings depending on if it was diagonal, down, right. This should make it O(mn).

He wanted to use a numpy array.

### Theoretical Analysis - Unrestricted Alignment

#### Time 

*Fill me in*

#### Space

*Fill me in*

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
