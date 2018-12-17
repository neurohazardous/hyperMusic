'''
This script generates the table.dat to run all analysis on the Graham cluster
'''
import random

pairs = ['P03', 'P04', 'P05', 'P08', 'P09', 'P11']
pairs_subs = {'P03': ('3A', '4B') , 'P04': ('5B', '6B'), 'P05': ('5A', '5B'),
              'P08':('8A', '8B'), 'P09': ('9A', '9B'), 'P11': ('11A', '11B')}
delays = ['10', '96', '480']
freqs = ['delta', 'theta', 'alpha', 'beta', 'gamma']
already_tested = []
for i in range(30):
    i_a = random.randint(0, 5)
    i_b = random.randint(0, 5)
    while i_a == i_b:
        i_b = random.randint(0, 5)
    pair_combination = str(i_a) + (str(i_b))

    while pair_combination in already_tested:
        i_a = random.randint(0, 5)
        i_b = random.randint(0, 5)
        while i_a == i_b:
            i_b = random.randint(0, 5)
        pair_combination = str(i_a) + (str(i_b))

    already_tested.append(str(i_a) + (str(i_b)))

    for d in delays:
        for freq_s in freqs:
            for freq_t in freqs:
                with open('5-STE_scrambled_pairs/table.dat', 'a') as the_file:
                    the_file.write('module load python/2.7.14; module load scipy-stack; python /home/horozco/5-STE_scrambled_pairs/5-analysis_scrambled_ste.py '
                           + pairs[i_a] + ' ' + pairs[i_b] + ' ' + freq_s + ' ' + freq_t + ' ' + pairs_subs[pairs[i_a]][0] + ' ' + pairs_subs[pairs[i_b]][1] + ' ' + d + '\n')

