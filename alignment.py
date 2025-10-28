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
    dist = {}
    prev = {}
    for i in range(-1, len(seq1)):
        dist[i] = {}
        for j in range(-1, len(seq2)):
            calculate_distance(seq1, seq2, indel_penalty, dist, prev, i, j, match_award, sub_penalty)
    edits = get_edit_sequence(seq1, seq2, prev)
    aligned_seq1, aligned_seq2 = get_sequence_alignment(seq1, seq2, gap, edits)
    return dist[len(seq1)-1][len(seq2)-1], aligned_seq1, aligned_seq2

def get_sequence_alignment(seq1, seq2, gap, edits):
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

def calculate_distance(seq1, seq2, indel_penalty, dist, prev, i, j, match_award, sub_penalty):
    if i == -1 and j == -1:
        dist[i][j] = 0
    elif i == -1:
        dist[i][j] = dist[i][j-1] + indel_penalty
        prev[i, j] = (i, j-1)
    elif j == -1:
        dist[i][j] = dist[i-1][j] + indel_penalty
        prev[i, j] = (i-1, j)
    else:
        up_scenario = dist[i-1][j] + indel_penalty
        left_scenario = dist[i][j-1] + indel_penalty
        diag_modifier = match_award if seq1[i] == seq2[j] else sub_penalty
        diag_scenario = dist[i-1][j-1] + diag_modifier
        min_scenario = min(up_scenario, left_scenario, diag_scenario)
        dist[i][j] = min_scenario
        if min_scenario == diag_scenario:
            prev[i, j] = (i-1, j-1)
        elif min_scenario == left_scenario:
            prev[i, j] = (i, j-1)
        else:
            prev[i, j] = (i-1, j)

