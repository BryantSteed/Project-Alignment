import math

def align(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap_open_penalty=0,
        gap='-',
) -> tuple[float, str | None, str | None]:
    """
        Align seq1 against seq2 using Needleman-Wunsch
        Put seq1 on left (j) and seq2 on top (i)
        => matrix[i][j]
        :param seq1: the first sequence to align; should be on the "left" of the matrix
        :param seq2: the second sequence to align; should be on the "top" of the matrix
        :param match_award: how many points to award a match
        :param indel_penalty: how many points to award a gap in either sequence
        :param sub_penalty: how many points to award a substitution
        :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
        :param gap_open_penalty: how much it costs to open a gap. If 0, there is no gap_open penalty
        :param gap: the character to use to represent gaps in the alignment strings
    """
    calculator = AlignmentCalculator(seq1, seq2, banded_width, 
                                     indel_penalty, match_award, 
                                     sub_penalty, gap_open_penalty)
    dist, prev = calculator.calculate_matrix()
    edits = get_edit_sequence(seq1, seq2, prev)
    aligned_seq1, aligned_seq2 = get_sequence_alignment(seq1, seq2, gap, edits)
    score = math.inf if aligned_seq1 is None else dist[len(seq1)-1][len(seq2)-1]
    return score, aligned_seq1, aligned_seq2

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
        self.dist = {i: {} for i in range(-1, len(self.seq1))}
        self.prev = {}
        self.mode = {i : {} for i in range(-1, len(self.seq1))}
    
    def calculate_matrix(self):
        for i, j in self.position_iterator():
            self.calculate_distance(i, j)
        return self.dist, self.prev

    def position_iterator(self):
        if self.banded_width == -1:
            for i in range(-1, len(self.seq1)):
                for j in range(-1, len(self.seq2)):
                    yield i, j
        else:
            for i in range(-1, len(self.seq1)):
                j_start = max(-1, i - self.banded_width)
                j_end = min(len(self.seq2) - 1, i + self.banded_width)
                for j in range(j_start, j_end + 1):
                    yield i, j
    
    def calculate_distance(self, i, j):
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
    
    def get_best_distance(self, i, j):
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
    
    def get_scenarios(self, i, j):
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

def get_sequence_alignment(seq1, seq2, gap, edits):
    if (not edits) and (len(seq1) != 0 and len(seq2) != 0):
        return None, None
    aligned_seq1, aligned_seq2 = "", ""
    seq1_shifts, seq2_shifts = 0, 0
    for edit in edits:
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
    edits = rev_edits[::-1]
    return edits