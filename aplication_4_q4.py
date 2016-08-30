import project_4 as pr
import application_4_data as data
import random
import matplotlib.pyplot as plt


# Question 4
print
print("Question 4 answer")

human_seq = data.read_protein(data.HUMAN_EYELESS_URL)
fly_seq = data.read_protein(data.FRUITFLY_EYELESS_URL)
scoring_matrix = data.read_scoring_matrix(data.PAM50_URL)

def generate_null_distrubution(seq_x, seq_y, scoring_matrix, num_trials):

    scoring_distribution = {}

    for dummy_i in xrange(num_trials):
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

number_of_trials = 1000
plot_data = generate_null_distrubution(human_seq, fly_seq, scoring_matrix, number_of_trials)
#print plot_data
x_plot_data = plot_data.keys()
#values_plot_data = plot_data.values()
y_plot_data = [value / float(number_of_trials) for value in plot_data.values()]
#print x_plot_data
#print values_plot_data
#print y_plot_data

plt.bar(x_plot_data, y_plot_data, 1, color='blue')
#plt.xscale('log')
#plt.yscale('log')
plt.grid(True)

plt.xlabel('score of local alignment')
plt.ylabel('fraction of total trials corsponding to score')
plt.title('normalizes scoring distribution of local alignment of sequence of Human Eyeless Protein and random permutation of Fruitfly Eyeless Protein')
plt.show()









