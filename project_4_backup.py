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
    """
    compute_alignment_matrix('A', 'A', {
    'A': {'A': 6, 'C': 2, '-': -4,
    'T': 2, 'G': 2}, 'C': {'A': 2,
    'C': 6, '-': -4, 'T': 2, 'G': 2},
    '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
    'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
    'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}
    }, False)
    expected [[0, 0], [0, 6]]
    """

    seq_x_len = len(seq_x)
    seq_y_len = len(seq_y)

    alignment_matrix = [[0]]

    for seq_x_letter in xrange(seq_x_len):
        previous_score = alignment_matrix[seq_x_letter][0]
        score = max(scoring_matrix[seq_x[seq_x_letter]]['-'], 0)
        if global_flag:
            score = scoring_matrix[seq_x[seq_x_letter]]['-']
        alignment_matrix.append([previous_score + score])

    for seq_y_letter in xrange(seq_y_len):
        print('y letter: ' + str(seq_y_letter))
        previous_score = alignment_matrix[0][seq_y_letter]
        print('prv score: ' + str(previous_score))
        score = max(scoring_matrix['-'][seq_y[seq_y_letter]], 0)
        if global_flag:
            score = scoring_matrix['-'][seq_y[seq_y_letter]]
        print('score: ' + str(score))
        alignment_matrix[0].append(previous_score + score)
    print
    for row in alignment_matrix:
        print(row)
    print
    print('nested loop')
    print
    for seq_x_letter in xrange(seq_x_len):
        for seq_y_letter in xrange(seq_y_len):

            x_y_score = scoring_matrix[seq_x[seq_x_letter]][seq_y[seq_y_letter]] + alignment_matrix[seq_x_letter][seq_y_letter]
            x_dash_score = scoring_matrix[seq_x[seq_x_letter]]['-'] + alignment_matrix[seq_x_letter][seq_y_letter +1]
            dash_y_score = scoring_matrix['-'][seq_y[seq_y_letter]] + alignment_matrix[seq_x_letter +1][seq_y_letter]

            score = max(x_y_score, x_dash_score, dash_y_score, 0)
                                          
            alignment_matrix[seq_x_letter +1].append(score)

    print
    for row in alignment_matrix:
        print(row)
    #return alignment_matrix

compute_alignment_matrix('cat', 'batcatdog', {
    '-': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    'a': {'-': -1, 'a': 2, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    'c': {'-': -1, 'a': -1, 'c': 2, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    'b': {'-': -1, 'a': -1, 'c': -1, 'b': 2, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    'd': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': 2, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    'g': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': 2, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    'o': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': 2, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
    't': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': 2, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}},
    False)

"""
expected [
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 2, 1, 0, 0, 0, 0],
[0, 0, 2, 1, 1, 4, 3, 2, 1, 0],
[0, 0, 1, 4, 3, 3, 6, 5, 4, 3]]

but received [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, -1, -1, 2, 1, 0, -1, -1, -1], [0, -1, 1, 0, 1, 4, 3, 2, 1, 0], [0, -1, 0, 3, 2, 3, 6, 5, 4, 3]]


compute_alignment_matrix('A', 'A', {
    'A': {'A': 6, 'C': 2, '-': -4,
    'T': 2, 'G': 2}, 'C': {'A': 2,
    'C': 6, '-': -4, 'T': 2, 'G': 2},
    '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
    'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
    'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}
    }, False)
"""
