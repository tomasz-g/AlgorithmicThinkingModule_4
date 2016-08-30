import project_4 as pr
import application_4_data as data


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
