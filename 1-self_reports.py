""" Hypermusic Self reports

This script intakes a CSV report data (specifically the MAQ questionnaire) and
compares the pre and post scores of how much the like each other and how much
they trust each other. Becuase of the nature of the questionnaire, we are
only including the last three questions for the pre and post:
4. I would like to become friends with my music partner
5. If my music partner needed help, I would help them
6. I would trust my music partner with a secret

To find this file, go to...
My drive > hyperMusic > hM_docs > hM_experiment > hM_questionnaires > 5.MAQ

Written by Hector D Orozco Perez

Warning: csv should be arranged by...
Pre_Participant A
Post_Participant B
Pre_Participant A
Post_Participant B

Last updated: 12/03/18
"""

import csv

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


# Get the csv file. Look at docstring for specific formatting!
with open('/Users/hectorOrozco/Desktop/hM_analysis/'
          'data/MAQ_responses120318.csv', 'rb') as f:
    reader = csv.reader(f)
    MAQ_raw = list(reader)

# Run the next loop for each pair of participants
participants = 3
part_id = []
pre_scores = np.empty(participants)
post_scores = np.empty(participants)

for p in range(participants):
    part_id.append('P0' + str(p+1))
    # Get the mean of each pre-score separately for each person
    PA_pre = np.mean(np.asarray(MAQ_raw[(4*p) + 1][6:9], dtype=float))
    PB_pre = np.mean(np.asarray(MAQ_raw[(4*p) + 2][6:9], dtype=float))

    # Get the mean of each post-score separately for each person
    PA_post = np.mean(np.asarray(MAQ_raw[(4*p) + 3][6:9], dtype=float))
    PB_post = np.mean(np.asarray(MAQ_raw[(4*p) + 4][6:9], dtype=float))

    #Get aggregated pre and post scores
    pre_scores[p] = np.mean([PA_pre, PB_pre])
    post_scores[p] = np.mean([PA_post, PB_post])

temp = np.concatenate((np.expand_dims(pre_scores, axis=1),
                         np.expand_dims(post_scores, axis=1)), axis=1)


# Create panda's data structure
scores = pd.DataFrame(data=temp,
                   index=part_id,
                   columns=['Pre', 'Post'])

# Plot them results
plt.clf()
sns.set_style("darkgrid")
ax = sns.pointplot(data=scores, color="#FF0000", ci = 'sd')
ax.plot(label = 'medium')
ax.set(ylabel='Scores (1-6)')
ax.figure.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                  'figs/1-self_reports.eps', format='eps', dpi=1000)





