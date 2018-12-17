'''
This script generates the table.dat to run all analysis on the Graham cluster
'''

pairs = ['P03', 'P04', 'P05', 'P08', 'P09', 'P11']
subjects = [('3A', '4B'), ('5B', '6B'), ('5A', '5B'), ('8A', '8B'), ('9A', '9B'), ('11A', '11B')]
delays = ['1', '5', '10', '48', '96', '144', '240', '480', '960']
freqs = ['delta', 'theta', 'alpha', 'beta', 'gamma']

for i in range(len(pairs)):
    for d in delays:
        with open('2-STE_graham/table.dat', 'a') as the_file:
            the_file.write(
                'module load python/2.7.14; module load scipy-stack; '
                'python /home/horozco/2-STE_graham/2-analysis_ste.py '
                + pairs[i] + ' ' + 'all' + ' ' + 'all' + ' ' + subjects[i][0] + ' ' +
                subjects[i][1] + ' ' + d + '\n')
        for freq_s in freqs:
            for freq_t in freqs:
                with open('2-STE_graham/table.dat', 'a') as the_file:
                    the_file.write('module load python/2.7.14; module load scipy-stack; python /home/horozco/2-STE_graham/2-analysis_ste.py '
                           + pairs[i] + ' ' + freq_s + ' ' + freq_t + ' ' + subjects[i][0] + ' ' + subjects[i][1] + ' ' + d + '\n')
