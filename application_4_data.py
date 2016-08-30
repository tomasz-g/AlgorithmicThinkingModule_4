"""
Provide code and solution for Application 4
"""

import math
import random
import urllib2
import matplotlib.pyplot as plt
#import alg_project4_solution as student

# data files
HUMAN_EYELESS_URL = 'C:/Users/tom/Desktop/Coursea/Project_4_data/HumanEyelessProtein.txt'
FRUITFLY_EYELESS_URL = 'C:/Users/tom/Desktop/Coursea/Project_4_data/FruitflyEyelessProtein.txt'
PAM50_URL = 'C:/Users/tom/Desktop/Coursea/Project_4_data/PAM50.txt'
CONSENSUS_PAX_URL = 'C:/Users/tom/Desktop/Coursea/Project_4_data/ConsensusPAXDomain.txt'
WORD_LIST_URL = 'C:/Users/tom/Desktop/Coursea/Project_4_data/words_list.txt'


###############################################
# provided code

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = open(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict




def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = open(filename)
    protein_seq = protein_file.read()
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = open(filename)
    
    # read in files as string
    words = word_file.read()
    
    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print "Loaded a dictionary with", len(word_list), "words"
    return word_list
