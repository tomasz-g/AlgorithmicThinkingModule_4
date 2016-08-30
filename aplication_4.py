import project_4 as pr
import application_4_data as data
import random


# Question 1
"""
Compute the local alignments of the sequences of HumanEyelessProtein and
FruitflyEyelessProtein using the PAM50 scoring matrix and enter the score
and local alignments for these two sequences.
"""
print
print("Question 1 answer")
human_seq = data.read_protein(data.HUMAN_EYELESS_URL)
fly_seq = data.read_protein(data.FRUITFLY_EYELESS_URL)
scoring_matrix = data.read_scoring_matrix(data.PAM50_URL)
human_fly_alignment_matrix_local = pr.compute_alignment_matrix(human_seq, fly_seq, scoring_matrix, False)

human_fly_local_alignment = pr.compute_local_alignment(human_seq, fly_seq, scoring_matrix, human_fly_alignment_matrix_local)
human_fly_local_alignment_score = human_fly_local_alignment[0]
human_local_alignment = human_fly_local_alignment[1]
fruitfly_local_alignment = human_fly_local_alignment[2]

print("Human and Fruitfly local alignment \n")
print("Length of Human alignment: " + str(len(human_local_alignment)))
print("Length of Fruitfly alignment: " + str(len(fruitfly_local_alignment)))
print("Score: " + str(human_fly_local_alignment_score))
print("Human     alignment: " + str(human_local_alignment))
print("Fruitfly  alignment: " + str(fruitfly_local_alignment))
print

# Question 2
"""
Load the file ConsensusPAXDomain. For each of the two sequences of the local
alignment computed in Question 1, do the following:

- Delete any dashes '-' present in the sequence.
- Compute the global alignment of this dash-less sequence with
  the ConsensusPAXDomain sequence.
- Compare corresponding elements of these two globally-aligned sequences
  (local vs. consensus) and compute the percentage of elements in these
  two sequences that agree.
"""
print
print("Question 2 answer")

Consensus_PAXDomain = data.read_protein(data.CONSENSUS_PAX_URL)

dash_less_human_local_alignment = ""
for letter in human_local_alignment:
    if letter != "-":
        dash_less_human_local_alignment += letter

dash_less_fly_local_alignment = ""
for letter in fruitfly_local_alignment:
    if letter != "-":
        dash_less_fly_local_alignment += letter

dash_less_human_local_alignment_PAX_seq__alignment_matrix_global = pr.compute_alignment_matrix(dash_less_human_local_alignment, Consensus_PAXDomain, scoring_matrix, True)
dash_less_human_local_alignment_PAX_seq__global_alignment = pr.compute_global_alignment(dash_less_human_local_alignment, Consensus_PAXDomain, scoring_matrix, dash_less_human_local_alignment_PAX_seq__alignment_matrix_global)

human_PAX_alignment_score = dash_less_human_local_alignment_PAX_seq__global_alignment[0]
human_global_alignment = dash_less_human_local_alignment_PAX_seq__global_alignment[1]
PAX_1_global_alignment = dash_less_human_local_alignment_PAX_seq__global_alignment[2]

human_global_alignment_length = len(human_global_alignment)
human_PAX_agree_percentage = float(0)
for letter in xrange(human_global_alignment_length):
    if human_global_alignment[letter] == PAX_1_global_alignment[letter]:
        human_PAX_agree_percentage += 1
human_PAX_agree_percentage = round((human_PAX_agree_percentage / human_global_alignment_length *100), 1)

dash_less_fly_local_alignment_PAX_seq__alignment_matrix_global = pr.compute_alignment_matrix(dash_less_fly_local_alignment, Consensus_PAXDomain, scoring_matrix, True)
dash_less_fly_local_alignment_PAX_seq__global_alignment = pr.compute_global_alignment(dash_less_fly_local_alignment, Consensus_PAXDomain, scoring_matrix, dash_less_fly_local_alignment_PAX_seq__alignment_matrix_global)

fly_PAX_alignment_score = dash_less_fly_local_alignment_PAX_seq__global_alignment[0]
fly_global_alignment = dash_less_fly_local_alignment_PAX_seq__global_alignment[1]
PAX_2_global_alignment = dash_less_fly_local_alignment_PAX_seq__global_alignment[2]

fly_global_alignment_length = len(fly_global_alignment)
fly_PAX_agree_percentage = float(0)
for letter in xrange(fly_global_alignment_length):
    if fly_global_alignment[letter] == PAX_2_global_alignment[letter]:
        fly_PAX_agree_percentage += 1
fly_PAX_agree_percentage = round((fly_PAX_agree_percentage / fly_global_alignment_length *100), 1)

print("Global alignment of dash-less Human local alignment and PAX Domain sequence \n")
print("Length of Human alignment: " + str(human_global_alignment_length))
print("Length of PAX Domain alignment: " + str(len(PAX_1_global_alignment)))
print("Score: " + str(human_PAX_alignment_score))
print("Percentage of elements that agree: " + str(human_PAX_agree_percentage) + "%")
print("Human alignment: " + str(human_global_alignment))
print("PAX   alignment: " + str(PAX_1_global_alignment))
print

print("Global alignment of dash-less Fruitfly local alignment and PAX Domain sequence \n")
print("Length of Fruitfly alignment: " + str(fly_global_alignment_length))
print("Length of PAX Domain alignment: " + str(len(PAX_2_global_alignment)))
print("Score: " + str(fly_PAX_alignment_score))
print("Percentage of elements that agree: " + str(fly_PAX_agree_percentage) + "%")
print("Fruitfly alignment: " + str(fly_global_alignment))
print("PAX      alignment: " + str(PAX_2_global_alignment))
print


# Question 4
print
print("Question 4 answer")

def generate_null_distrubution(seq_x, seq_y, scoring_matrix, num_trials):

    scoring_distribution = {}

    for dummy_i in xrange(num_trils):
        rand_y = list(seq_y)
        random.shuffle(rand_y)
        rand_y = ''.join(rand_y)

        seq_x_rand_y_alignment_matrix = pr.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        seq_x_rand_y_local_alignment = pr.compute_local_alignment(seq_x, rand_y, scoring_matrix, seq_x_rand_y_alignment_matrix)
        seq_x_rand_y_maximum_score = seq_x_rand_y_local_alignment[0]
        if seq_x_rand_y_maximum_score in scoring_distribution.keys():
            scoring_distribution[seq_x_rand_y_maximum_score] += 1
        else:
            scoring_distribution[seq_x_rand_y_maximum_score] = 1
            
    return scoring_distribution








    
