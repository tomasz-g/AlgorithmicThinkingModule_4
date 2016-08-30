"""
Project 4

Two Matrix functions - return matrices
Two Alignment functions - return global and local alignments
"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """ (set([str]), int, int, int) -> dict

     Return a scoring matrix as a dictionary of dictionaries
     
     >>> build_scoring_matrix(set(['A', 'C', 'T', 'G']), 6, 2, -4)
     {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
     'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
     '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
     'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
     'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}     
     """

    scoring_matrix = {}
    scoring_matrix['-'] = {}
    scoring_matrix['-']['-'] = dash_score
    for key in alphabet:
        scoring_matrix[key] = {}
        scoring_matrix[key]['-'] = dash_score
        scoring_matrix['-'][key] = dash_score
        for value in alphabet:
            score = off_diag_score
            if key == value:
                score = diag_score
            scoring_matrix[key][value] = score
    return scoring_matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """ (str, str, dict, bool) -> list

    Return a alignment matrix as a list of lists

    >>> compute_alignment_matrix('A', 'A',
    {'A': {'A': 6, 'C': 2, '-': -4,
    'T': 2, 'G': 2}, 'C': {'A': 2,
    'C': 6, '-': -4, 'T': 2, 'G': 2},
    '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
    'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
    'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}},
    False)
    [[0, 0], [0, 6]]

    >>> compute_alignment_matrix('A', 'A',
    {'A': {'A': 6, 'C': 2, '-': -4,
    'T': 2, 'G': 2}, 'C': {'A': 2,
    'C': 6, '-': -4, 'T': 2, 'G': 2},
    '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
    'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
    'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}},
    True)
    [[0, -4], [-4, 6]]
    """
    
    seq_x_len = len(seq_x)
    seq_y_len = len(seq_y)
    
    # if 'global_flag' = False, all negative values computed for 'alignment_matrix' are changes to 0 
    flag = 0
    if global_flag:
        flag = float('-inf')
        
    alignment_matrix = [[0]]

    for seq_x_letter in xrange(seq_x_len):
        previous_score = alignment_matrix[seq_x_letter][0]
        score = max(scoring_matrix[seq_x[seq_x_letter]]['-'], flag)
        alignment_matrix.append([previous_score + score])

    for seq_y_letter in xrange(seq_y_len):
        previous_score = alignment_matrix[0][seq_y_letter]
        score = max(scoring_matrix['-'][seq_y[seq_y_letter]], flag)
        alignment_matrix[0].append(previous_score + score)
    
    for seq_x_letter in xrange(seq_x_len):
        for seq_y_letter in xrange(seq_y_len):
            x_y_score = scoring_matrix[seq_x[seq_x_letter]][seq_y[seq_y_letter]] + alignment_matrix[seq_x_letter][seq_y_letter]
            x_dash_score = scoring_matrix[seq_x[seq_x_letter]]['-'] + alignment_matrix[seq_x_letter][seq_y_letter +1]
            dash_y_score = scoring_matrix['-'][seq_y[seq_y_letter]] + alignment_matrix[seq_x_letter +1][seq_y_letter]
            score = max(x_y_score, x_dash_score, dash_y_score, flag)                   
            alignment_matrix[seq_x_letter +1].append(score)

    return alignment_matrix

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """ (str, str, dict, list) -> tuple

    Return score of
    """
    
    seq_x_len = len(seq_x)
    seq_y_len = len(seq_y)
    align_x = ""
    align_y = ""
    score = alignment_matrix[-1][-1] 

    while seq_x_len > 0 and seq_y_len > 0:
        if alignment_matrix[seq_x_len][seq_y_len] == alignment_matrix[seq_x_len -1][seq_y_len -1] + scoring_matrix[seq_x[seq_x_len -1]][seq_y[seq_y_len -1]]:
            align_x = seq_x[seq_x_len -1] + align_x
            align_y = seq_y[seq_y_len -1] + align_y
            seq_x_len -= 1
            seq_y_len -= 1

        else:            
            if alignment_matrix[seq_x_len][seq_y_len] == alignment_matrix[seq_x_len -1][seq_y_len] + scoring_matrix[seq_x[seq_x_len -1]]["-"]:
                align_x = seq_x[seq_x_len -1] + align_x
                align_y = "-" + align_y
                seq_x_len -= 1
            else:
                align_x = "-" + align_x
                align_y = seq_y[seq_y_len -1] + align_y
                seq_y_len -= 1

    while seq_x_len > 0:
        align_x = seq_x[seq_x_len -1] + align_x
        align_y = "-" + align_y
        seq_x_len -= 1

    while seq_y_len > 0:
        align_x = "-" + align_x
        align_y = seq_y[seq_y_len -1] + align_y
        seq_y_len -= 1

    return (score, align_x, align_y)

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """ (str, str, dict, list) -> tuple

    Return score of
    """
    seq_x_idx = 0
    seq_y_idx = 0
    align_x = ""
    align_y = ""

    score = 0
    for row in xrange(len(alignment_matrix)):
        for col in xrange(len(alignment_matrix[row])):
            matrix_score = alignment_matrix[row][col]
            if matrix_score >= score:
                score = matrix_score
                seq_x_idx = row
                seq_y_idx = col

    while alignment_matrix[seq_x_idx][seq_y_idx] != 0:
        if alignment_matrix[seq_x_idx][seq_y_idx] == alignment_matrix[seq_x_idx -1][seq_y_idx -1] + scoring_matrix[seq_x[seq_x_idx -1]][seq_y[seq_y_idx -1]]:
            align_x = seq_x[seq_x_idx -1] + align_x
            align_y = seq_y[seq_y_idx -1] + align_y
            seq_x_idx -= 1
            seq_y_idx -= 1

        else:            
            if alignment_matrix[seq_x_idx][seq_y_idx] == alignment_matrix[seq_x_idx -1][seq_y_idx] + scoring_matrix[seq_x[seq_x_idx -1]]["-"]:
                align_x = seq_x[seq_x_idx -1] + align_x
                align_y = "-" + align_y
                seq_x_idx -= 1
            else:
                align_x = "-" + align_x
                align_y = seq_y[seq_y_idx -1] + align_y
                seq_y_idx -= 1

    return (score, align_x, align_y)




