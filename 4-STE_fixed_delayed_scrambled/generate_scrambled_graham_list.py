'''
This script generates the table.dat to run all analysis on the Graham cluster
'''
import random

pairs = ['P03', 'P04', 'P05', 'P08', 'P09', 'P11']
pairs_subs = {'P03': ('3A', '4B') , 'P04': ('5B', '6B'), 'P05': ('5A', '5B'),
              'P08':('8A', '8B'), 'P09': ('9A', '9B'), 'P11': ('11A', '11B')}
delays = ['3', '30', '150']
freqs = ['delta', 'theta', 'alpha', 'beta', 'gamma']

for i_pairs in range(1, 2):
    for d in delays:
        for freq_s in freqs:
            for freq_t in freqs:
                with open('8-STE_fixed_delayed_scrambled/table.dat', 'a') as the_file:
                    if i_pairs == 5:
                        i_second = 0
                    elif i_pairs == 1:
                        i_second = 0
                    else:
                        i_second = i_pairs + 1
                    the_file.write('module load python/2.7.14; module load scipy-stack; python /home/horozco/8-STE_fixed_delayed_scrambled/8-analysis_scrambled_ste.py '
                           + pairs[i_pairs] + ' ' + pairs[i_second] + ' ' + freq_s + ' ' + freq_t + ' ' + pairs_subs[pairs[i_pairs]][0] + ' ' + pairs_subs[pairs[i_second]][1] + ' ' + d + '\n')