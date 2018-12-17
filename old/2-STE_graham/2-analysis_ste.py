'''
== License ===================================================================
This file is part of the project hyperMusic. All of
hyperMusic code is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
hyperMusic is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with hyperMusic.
If not, see <https://www.gnu.org/licenses/>.

== Remarks ===================================================================
This module includes all the necessary functions to calculate symbolic
transfer entropy using NumPy and Matplotlib for the project
hyperMusic.

This code was translated from a MATLAB implementation of this algorithm done
by the BIAPT lab (https://github.com/BIAPT/EEGapp).

See "Lee, U., Blain-Moraes, S., & Mashour, G. A. (2015).
Assessing levels of consciousness with symbolic analysis.
Phil. Trans. R. Soc. A, 373(2034), 20140117." for a nice intro to the topic.

For a more detailed description of the algorithm, see...
(1) Schreiber, T. (2000). Measuring information transfer.
Physical review letters, 85(2), 461.

(2) Staniek, M., & Lehnertz, K. (2008). Symbolic transfer entropy.
Physical Review Letters, 100(15), 158101.

Python implementation: Hector D Orozco Perez

This code is meant to run on the Graham (Compute Canada) system using in a
serial farm

== hyperMusic: STE Analysis ==================================================
1) Get parameters from terminal input
2) Get data from .mat files
3) Baseline time-frequency decomposition (for normalization)
4) Trial time-frequency decomposition
5) Symbolic Transfer Entropy computation
'''

import sys
import time

import numpy as np
import scipy as sp
from scipy import signal
import scipy.io as sio
import pandas as pd


def symb_transfer_entropy(source, target, delta, k):
    '''
    Calculates symbolic transfer entropy of source to target. The
    calculation
    is done column (the whole "set" of possible events) to row (specific
    event)
    Args:
        source (numpy.ndarray): 2D array the form of [time, 1]
        target (numpy.ndarray): 2D array the form of [time, 1]
        delta (int): the lag between the signal and the prediction
        k (int): embedding dimension for the symbolization process
    Returns:
        symb_trans_entropy (float): raw symbolic transfer entropy values
    '''

    # Wrangle data into [time,k] matrix
    N = len(source)
    delayed_data_s = np.array(
            [source[i:i + k] for i in range(N - (k - 1))])
    delayed_data_t = np.array(
            [target[i:i + k] for i in range(N - (k - 1))])

    # Symbolize the data by ranking each of the row's columns and
    # concatenating those numbers using 10 and its powers
    symbols_s = np.squeeze(np.argsort(delayed_data_s, axis=1) + 1)
    integer_s = np.array(
            [symbols_s[:, i] * 10 ** (k - 1 - i) for i in range(k)])
    integer_s = integer_s.T.sum(axis=1)[:, None]

    symbols_t = np.squeeze(np.argsort(delayed_data_t, axis=1) + 1)
    integer_t = np.array(
            [symbols_t[:, i] * 10 ** (k - 1 - i) for i in range(k)])
    integer_t = integer_t.T.sum(axis=1)[:, None]

    # Calculate probabilities needed for STE
    integer_both = np.concatenate((integer_t, integer_s),
                                  axis=1)  # target, source
    int_future = integer_both[delta:, :]
    int_past = integer_both[0:-delta, :]

    # We first define an internal function to estimate probability
    def estimate_prob(data):
        '''
        Estimate probability given the data provided (column is whole
         universe of probabilities, row is specific event)
        Args:
            data (numpy.ndarray): 1D array with all int symbols
        Returns:
            prob (numpy.ndarray): probability of each event
        '''
        count = np.zeros((data.shape[0], 1))
        b = np.unique(data)  #discrete/non-overlapping events

        def index_events(b_i):
            index = np.where(b_i == data)  #get the indexes of events
            count[index] = len(index[0])
            return count / float(data.shape[0])

        prob = np.squeeze(np.array(map(index_events, b)))[-1, :]
        return prob[:, None]

    # {Target future, Target past, Source past}
    int1 = int_future[:, 0] * 10 ** (k * 2) + int_past[:, 0] * 10 ** (
        k) + int_past[:, 1]

    # {Target past, Source past}
    int2 = int_past[:, 0] * 10 ** (k) + int_past[:, 1]

    # {Target future, Target past}
    int3 = int_future[:, 0] * 10 ** (k) + int_past[:, 0]

    # {Target past}
    int4 = int_past[:, 0]

    # Estimate actual probabilities
    p1 = estimate_prob(int1)  # P[target(n+delta), target(n), source(n)]
    p2 = estimate_prob(int2)  # P[target(n), source(n)]
    p3 = estimate_prob(int3)  # P[target(n+delta), target(n)]
    p4 = estimate_prob(int4)  # P[target(n)]

    # int1 gives us the set of unique possible outcomes
    _, unique_index = np.unique(int1, return_index=True)
    p1 = p1[unique_index, :]
    p2 = p2[unique_index, :]
    p3 = p3[unique_index, :]
    p4 = p4[unique_index, :]

    # Calculate symbolic transfer entropy
    symb_trans_entropy = (p1 * (np.log2(p1 * p4) - np.log2(p2 * p3))).sum()
    return symb_trans_entropy

# A little test to see it works as expected

import matplotlib.pyplot as plt
x = np.random.rand(500, 1)
y = np.zeros((500,1))
y[10:] = x[:-10]
plt.plot(x[0:20])
plt.plot(y[0:20])
plt.show()

delta = np.linspace(1, 20, 20).astype(int)
delta_plot_x2y = np.zeros(len(delta))
delta_plot_y2x = np.zeros(len(delta))

for i_d in range(len(delta)):
    delta_plot_x2y[i_d] = symb_transfer_entropy(x, y, delta[i_d], 3)
    delta_plot_y2x[i_d] = symb_transfer_entropy(y, x, delta[i_d], 3)

plt.plot(delta, delta_plot_x2y, label = 'source to target')
plt.plot(delta, delta_plot_y2x, label = 'target to source')
plt.legend()
plt.show()


def main():
    # == 1) Parameters =======================================================
    pair = str(sys.argv[1])
    freq_band_s = str(sys.argv[2])
    freq_band_t = str(sys.argv[3])
    sub_a = str(sys.argv[4])
    sub_b = str(sys.argv[5])
    delay = int(sys.argv[6])

    fsample = 2400.
    patches_to_remove = ['Basal Ganglia Left', 'Basal Ganglia Right',
                         'Limbic Left', 'Limbic Right','Cerebellum Left',
                         'Cerebellum Right', 'Cerebellum Mid']

    freq_all = {'delta': [1., 3.], 'theta': [4., 7.], 'alpha': [8., 12.],
                   'beta': [13., 29.], 'gamma': [30., 45.], 'all': False}
    freq_bound_s = freq_all[freq_band_s]
    freq_bound_t = freq_all[freq_band_t]

    # == 2) Get data from .mat files =========================================
    # Sub a -> arbitrary who it is
    # If ran in Graham, change directory to /home/horozco
    temp = sio.loadmat('/home/horozco/2-STE_graham/' + pair + '/' + sub_a +
                       '_trials_labels.mat')
    sub_a_trials_labels = [temp['trials_labels'][n][0][0].encode("utf-8")
                           for n in range(21)]
    del temp

    temp = sio.loadmat('/home/horozco/2-STE_graham/' + pair + '/' + sub_a +
                       '_source_labels.mat')
    sub_a_patch_labels = [temp['source_labels'][n][0][0].encode("utf-8")
                          for n in range(19)]
    [sub_a_patch_labels.remove(n) for n in patches_to_remove]
    patch_labels = sub_a_patch_labels
    del temp

    temp = sio.loadmat('/home/horozco/2-STE_graham/' + pair + '/' + sub_a +
                       '_sources.mat')
    sub_a_trials_data = [temp['sources'][n][0] for n in range(21)]
    del temp

    # Sub b -> arbitrary who it is
    temp = sio.loadmat('/home/horozco/2-STE_graham/' + pair + '/' + sub_b +
                       '_trials_labels.mat')
    sub_b_trials_labels = [temp['trials_labels'][n][0][0].encode("utf-8")
                           for n in range(21)]
    del temp

    temp = sio.loadmat('/home/horozco/2-STE_graham/' + pair + '/' + sub_b +
                       '_sources.mat')
    sub_b_trials_data = [temp['sources'][n][0] for n in range(21)]
    del temp

    # == 3) Baseline TF decomposition (for normalization) ====================
    # == Baseline: source frequency band =====================================
    # Get transition badnwdith using MNE's rule of thumb
    if freq_bound_s: # all is only compared to itself
        freq_l_trans_s = min(max(freq_bound_s[0] * 0.25, 2), freq_bound_s[0])
        freq_h_trans_s = min(max(freq_bound_s[1] * 0.25, 2.),
                             fsample / 2. - freq_bound_s[1])
        # Compute the filter
        filter_nyq = fsample / 2.0
        filter_s_width = np.mean(np.array([freq_l_trans_s, freq_h_trans_s]))
        filter_s_n, _ = signal.kaiserord(30.0, filter_s_width / filter_nyq)
        filter_s_taps = signal.firwin(filter_s_n,
                                      np.array(freq_bound_s),
                                      width=filter_s_width, window='blackman',
                                      pass_zero=False, fs=fsample)

    # == Sub a ===============================================================
    # Get rid of the deep brain structures
    sub_a_baseline_sfreq = sub_a_trials_data[-1]
    sub_a_baseline_sfreq = sub_a_baseline_sfreq[
                           [0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13], :]

    # Apply filter (forward and backward  using filtfilt)
    if freq_bound_s:
        sub_a_base_s_filt = signal.filtfilt(filter_s_taps, 1,
                                            sub_a_baseline_sfreq,
                                            axis = 1, padlen = None)
    else:
        sub_a_base_s_filt = sub_a_baseline_sfreq # Not filtered; broadband EEG

    # Get power from Hilbert Transform
    orig_n = sub_a_base_s_filt.shape[1]
    hilt_n = sp.fftpack.next_fast_len(orig_n) # Make HT go faster
    sub_a_base_s_env = signal.hilbert(sub_a_base_s_filt, axis = 1,
                                      N = hilt_n)
    sub_a_base_s_env = sub_a_base_s_env[:, :orig_n] # Crop out the padding
    sub_a_base_s_env_power = np.abs(sub_a_base_s_env)**2

    # Get rid of padding (3s at the begining and at the end) to take care of
    # HT and filtering edge artifacts
    padding = int(3 * fsample)
    sub_a_base_s_env_power = sub_a_base_s_env_power[:, padding:-padding]
    sub_a_baseline_s = sub_a_base_s_env_power.mean(axis = 1)

    # == Sub b ===============================================================
    # Get rid of the deep brain structures
    sub_b_baseline_sfreq = sub_b_trials_data[-1]
    sub_b_baseline_sfreq = sub_b_baseline_sfreq[
                           [0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13], :]

    # Apply filter (forward and backward  using filtfilt)
    if freq_bound_t:
        sub_b_base_s_filt = signal.filtfilt(filter_s_taps, 1,
                                            sub_b_baseline_sfreq,
                                            axis=1, padlen=None)
    else:
        sub_b_base_s_filt = sub_b_baseline_sfreq # Not filtered; broadband EEG

    # Get power from Hilbert Transform
    orig_n = sub_b_base_s_filt.shape[1]
    hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
    sub_b_base_s_env = signal.hilbert(sub_b_base_s_filt, axis=1,
                                      N=hilt_n)
    sub_b_base_s_env = sub_b_base_s_env[:, :orig_n]  # Crop out the padding
    sub_b_base_s_env_power = np.abs(sub_b_base_s_env) ** 2

    # Get rid of padding (3s at the begining and at the end) to take care of
    # HT and filtering edge artifacts
    padding = int(3 * fsample)
    sub_b_base_s_env_power = sub_b_base_s_env_power[:, padding:-padding]
    sub_b_baseline_s = sub_b_base_s_env_power.mean(axis=1)

    if freq_band_s == freq_band_t:
        sub_a_baseline_t = sub_a_baseline_s
        sub_b_baseline_t = sub_b_baseline_s
    else:
        # == Baseline: target frequency band =================================
        # Get transition bandwidth using MNE's rule of thumb
        freq_l_trans_t = min(max(freq_bound_t[0] * 0.25, 2), freq_bound_t[0])
        freq_h_trans_t = min(max(freq_bound_t[1] * 0.25, 2.),
                             fsample / 2. - freq_bound_t[1])
        # Compute the filter
        filter_nyq = fsample / 2.0
        filter_t_width = np.mean(np.array([freq_l_trans_t, freq_h_trans_t]))
        filter_t_n, _ = signal.kaiserord(30.0, filter_t_width / filter_nyq)
        filter_t_taps = signal.firwin(filter_t_n,
                                      np.array(freq_bound_t),
                                      width=filter_t_width, window='blackman',
                                      pass_zero=False, fs=fsample)

        # == Sub a ===========================================================
        # Get rid of the deep brain structures
        sub_a_baseline_tfreq = sub_a_trials_data[-1]
        sub_a_baseline_tfreq = sub_a_baseline_tfreq[
                               [0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13], :]

        # Apply filter (forward and backward  using filtfilt)
        sub_a_base_t_filt = signal.filtfilt(filter_t_taps, 1,
                                            sub_a_baseline_tfreq,
                                            axis=1, padlen=None)

        # Get power from Hilbert Transform
        orig_n = sub_a_base_t_filt.shape[1]
        hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
        sub_a_base_t_env = signal.hilbert(sub_a_base_t_filt, axis=1,
                                          N=hilt_n)
        sub_a_base_t_env = sub_a_base_t_env[:, :orig_n]  # Crop out the padding
        sub_a_base_t_env_power = np.abs(sub_a_base_t_env) ** 2

        # Get rid of padding (3s at the begining and at the end) to take
        # care of HT and filtering edge artifacts
        padding = int(3 * fsample)
        sub_a_base_t_env_power = sub_a_base_t_env_power[:, padding:-padding]
        sub_a_baseline_t = sub_a_base_t_env_power.mean(axis=1)

        # == Sub b ===========================================================
        # Get rid of the deep brain structures
        sub_b_baseline_tfreq = sub_b_trials_data[-1]
        sub_b_baseline_tfreq = sub_b_baseline_tfreq[
                               [0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13], :]

        # Apply filter (forward and backward  using filtfilt)
        sub_b_base_t_filt = signal.filtfilt(filter_t_taps, 1,
                                            sub_b_baseline_tfreq,
                                            axis=1, padlen=None)

        # Get power from Hilbert Transform
        orig_n = sub_b_base_t_filt.shape[1]
        hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
        sub_b_base_t_env = signal.hilbert(sub_b_base_t_filt, axis=1,
                                          N=hilt_n)
        sub_b_base_t_env = sub_b_base_t_env[:, :orig_n]  # Crop out the padding
        sub_b_base_t_env_power = np.abs(sub_b_base_t_env) ** 2

        # Get rid of padding (3s at the begining and at the end) to take
        # care of HT and filtering edge artifacts
        padding = int(3 * fsample)
        sub_b_base_t_env_power = sub_b_base_t_env_power[:, padding:-padding]
        sub_b_baseline_t = sub_b_base_t_env_power.mean(axis=1)

    # == 4) Trial Time Frequency Decomposition ===============================
    for n_trial in range(len(sub_a_trials_labels)-1):
        start_time = time.time()
        print('Now processing trial %i' % (n_trial+1))
        # == Source frequency band ===========================================
        # == Sub a ===========================================================
        # Get rid of deep brain structures
        sub_a_trial = sub_a_trials_data[n_trial]
        sub_a_trial = sub_a_trial[[0, 1, 2, 3, 6, 7, 8, 9 , 10, 11, 12, 13], :]

        # Apply filter (forward and backward  using filtfilt)
        if freq_bound_s:
            sub_a_trial_s_filt = signal.filtfilt(filter_s_taps, 1,
                                                 sub_a_trial, axis=1,
                                                 padlen=None)
        else:
            sub_a_trial_s_filt = sub_a_trial # Not filtered; broadband EEG

        # Get power from Hilbert Transform
        orig_n = sub_a_trial_s_filt.shape[1]
        hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
        sub_a_trial_s_env = signal.hilbert(sub_a_trial_s_filt, axis=1,
                                          N=hilt_n)
        sub_a_trial_s_env = sub_a_trial_s_env[:, :orig_n]  # Crop out the padding
        sub_a_trial_s_env_power = np.abs(sub_a_trial_s_env) ** 2

        # Do baseline normalization (% change from baseline)
        sub_a_baseline_s_tiled = np.tile(sub_a_baseline_s,
                                         (sub_a_trial_s_env_power.shape[1], 1)).T
        sub_a_trial_s_norm = ((sub_a_trial_s_env_power-sub_a_baseline_s_tiled)
                               /sub_a_baseline_s_tiled)*100

        # Get rid of padding (3s at the begining and at the end) to take
        # care of HT and filtering edge artifacts
        padding = int(3 * fsample)
        sub_a_trial_s_norm = sub_a_trial_s_norm[:, padding:-padding]

        # Downsample HT
        sub_a_trial_s = signal.decimate(sub_a_trial_s_norm, 5, axis=1, zero_phase=True)

        # == Sub b ===========================================================
        # Get rid of deep brain structures
        sub_b_trial = sub_b_trials_data[n_trial]
        sub_b_trial = sub_b_trial[[0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13], :]

        # Apply filter (forward and backward  using filtfilt)
        if freq_bound_s:
            sub_b_trial_s_filt = signal.filtfilt(filter_s_taps, 1,
                                             sub_b_trial, axis=1, padlen=None)
        else:
            sub_b_trial_s_filt = sub_b_trial # Not filtered; broadband EEG

        # Get power from Hilbert Transform
        orig_n = sub_b_trial_s_filt.shape[1]
        hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
        sub_b_trial_s_env = signal.hilbert(sub_b_trial_s_filt, axis=1,
                                           N=hilt_n)
        sub_b_trial_s_env = sub_b_trial_s_env[:, :orig_n]
        sub_b_trial_s_env_power = np.abs(sub_b_trial_s_env)** 2

        # Do baseline normalization (% change from baseline)
        sub_b_baseline_s_tiled = np.tile(sub_b_baseline_s, (sub_b_trial_s_env_power.shape[1],1)).T
        sub_b_trial_s_norm = ((sub_b_trial_s_env_power-sub_b_baseline_s_tiled)
                              /sub_b_baseline_s_tiled) * 100

        # Get rid of padding (3s at the begining and at the end) to take
        # care of HT and filtering edge artifacts
        padding = int(3 * fsample)
        sub_b_trial_s_norm = sub_b_trial_s_norm[:, padding:-padding]

        # Downsample HT
        sub_b_trial_s = signal.decimate(sub_b_trial_s_norm, 5, axis=1,zero_phase=True)

        # == Target frequency band ===========================================
        if freq_band_s == freq_band_t:
            sub_a_trial_t = sub_a_trial_s
            sub_b_trial_t = sub_b_trial_s
        else:
            # == Sub a ===========================================================
            # Apply filter (forward and backward  using filtfilt)
            sub_a_trial_t_filt = signal.filtfilt(filter_t_taps, 1,
                                                 sub_a_trial, axis=1, padlen=None)

            # Get power from Hilbert Transform
            orig_n = sub_a_trial_t_filt.shape[1]
            hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
            sub_a_trial_t_env = signal.hilbert(sub_a_trial_t_filt, axis=1,
                                               N=hilt_n)
            sub_a_trial_t_env = sub_a_trial_t_env[:, :orig_n]
            sub_a_trial_t_env_power = np.abs(sub_a_trial_t_env) ** 2

            # Do baseline normalization (% change from baseline)
            sub_a_baseline_t_tiled = np.tile(sub_a_baseline_t,(sub_a_trial_t_env_power.shape[1],1)).T
            sub_a_trial_t_norm = ((sub_a_trial_t_env_power -sub_a_baseline_t_tiled)
                                  / sub_a_baseline_t_tiled) * 100

            # Get rid of padding (3s at the begining and at the end) to take
            # care of HT and filtering edge artifacts
            padding = int(3 * fsample)
            sub_a_trial_t_norm = sub_a_trial_t_norm[:, padding:-padding]

            # Downsample HT
            sub_a_trial_t = signal.decimate(sub_a_trial_t_norm, 5, axis=1,
                                            zero_phase=True)

            # == Sub b ===========================================================
            # Apply filter (forward and backward  using filtfilt)
            sub_b_trial_t_filt = signal.filtfilt(filter_t_taps, 1,
                                                 sub_b_trial, axis=1, padlen=None)

            # Get power from Hilbert Transform
            orig_n = sub_b_trial_t_filt.shape[1]
            hilt_n = sp.fftpack.next_fast_len(orig_n)  # Make HT go faster
            sub_b_trial_t_env = signal.hilbert(sub_b_trial_t_filt, axis=1,
                                               N=hilt_n)
            sub_b_trial_t_env = sub_b_trial_t_env[:, :orig_n]
            sub_b_trial_t_env_power = np.abs(sub_b_trial_t_env) ** 2

            # Do baseline normalization (% change from baseline)
            sub_b_baseline_t_tiled = np.tile(sub_b_baseline_t, (sub_b_trial_t_env_power.shape[1], 1)).T
            sub_b_trial_t_norm = ((sub_b_trial_t_env_power-sub_b_baseline_t_tiled)
                                  / sub_b_baseline_t_tiled) * 100

            # Get rid of padding (3s at the begining and at the end) to take
            # care of HT and filtering edge artifacts
            padding = int(3 * fsample)
            sub_b_trial_t_norm = sub_b_trial_t_norm[:, padding:-padding]

            # Downsample HT
            sub_b_trial_t = signal.decimate(sub_b_trial_t_norm, 5, axis=1,
                                            zero_phase=True)

        # == 5) Symbolic Transfer Entropy ====================================
        temp_columns = ['sub_source', 'sub_target',
                                         'sub_leader', 'sub_follower',
                                         'freq_source', 'freq_target',
                                         'neuro_source', 'neuro_target',
                                         'duo_type', 'piece', 'trial',
                                         'pair', 'ste']
        temp_data = np.zeros(((sub_a_trial_s.shape[0] + sub_b_trial_s.shape[0]) *
                              (sub_a_trial_t.shape[0] + sub_b_trial_t.shape[0]),
                              len(temp_columns)))
        temp_df = pd.DataFrame(data=temp_data,
                                columns=temp_columns)
        # Important distinction: sub_follower and sub_leader refer to who
        # was piano I during the trial.
        # They remain constant at the level of the trial. sub_source and
        # sub_target refer to the direction
        # of information flow being analyze (leader to leader, follower to leader)


        temp_df['freq_source'] = freq_band_s
        temp_df['freq_target'] = freq_band_t
        temp_df['pair'] = pair
        temp_df['trial'] = sub_a_trials_labels[n_trial].lower()[-1]
        if int(pair[-2:]) > 9:
            temp_df['duo_type'] = sub_a_trials_labels[n_trial].lower()[4]
            temp_df['piece'] = sub_a_trials_labels[n_trial].lower()[5:8]
        else:
            temp_df['duo_type'] = sub_a_trials_labels[n_trial].lower()[3]
            temp_df['piece'] = sub_a_trials_labels[n_trial].lower()[4:7]

        trial_label = [sub_a_trials_labels[n_trial].lower(),
                       sub_b_trials_labels[n_trial].lower()]
        if trial_label[0][7] == 'a' and int(pair[-2:]) <= 9:
            leader_trial_s = sub_a_trial_s
            leader_trial_t = sub_a_trial_t
            follower_trial_s = sub_b_trial_s
            follower_trial_t = sub_b_trial_t
            temp_df['sub_leader'] = sub_a
            temp_df['sub_follower'] = sub_b
            leader = sub_a
            follower = sub_b
        elif trial_label[0][8] == 'a' and int(pair[-2:]) > 9:
            leader_trial_s = sub_a_trial_s
            leader_trial_t = sub_a_trial_t
            follower_trial_s = sub_b_trial_s
            follower_trial_t = sub_b_trial_t
            temp_df['sub_leader'] = sub_a
            temp_df['sub_follower'] = sub_b
            leader = sub_a
            follower = sub_b
        elif trial_label[1][7] == 'b' and int(pair[-2:]) <= 9:
            leader_trial_s = sub_b_trial_s
            leader_trial_t = sub_b_trial_t
            follower_trial_s = sub_a_trial_s
            follower_trial_t = sub_a_trial_t
            temp_df['sub_follower'] = sub_a
            temp_df['sub_leader'] = sub_b
            leader = sub_b
            follower = sub_a
        elif trial_label[1][8] == 'b' and int(pair[-2:]) > 9:
            leader_trial_s = sub_b_trial_s
            leader_trial_t = sub_b_trial_t
            follower_trial_s = sub_a_trial_s
            follower_trial_t = sub_a_trial_t
            temp_df['sub_follower'] = sub_a
            temp_df['sub_leader'] = sub_b
            leader = sub_b
            follower = sub_a
        else:
            print 'Unable to determine leader and follower, check trial %s' % trial_label
            break

        for n_source in range(len(patch_labels)):
            print("--- Starting %s ---" % patch_labels[n_source])
            temp_df.loc[(48*n_source):(47*(n_source+1)+n_source), 'neuro_source'] = patch_labels[n_source]
            temp_df.loc[(48*n_source):(47*(n_source+1)+n_source), 'neuro_target'] = (patch_labels * 4)

            # Leader to leader
            temp_df.loc[(48*n_source):(48*n_source+11), 'ste'] = [symb_transfer_entropy(leader_trial_s[n_source, :], leader_trial_t[i, :], delay, 3) for i in range(leader_trial_t.shape[0])]
            temp_df.loc[(48*n_source):(48*n_source+11), 'sub_source'] = leader
            temp_df.loc[(48*n_source):(48*n_source+11), 'sub_target'] = leader

            # Leader to follower
            temp_df.loc[(48*n_source+12):(48*n_source+23), 'ste'] = [symb_transfer_entropy(leader_trial_s[n_source, :], follower_trial_t[i, :], delay, 3) for i in range(leader_trial_t.shape[0])]
            temp_df.loc[(48*n_source+12):(48*n_source+23), 'sub_source'] = leader
            temp_df.loc[(48*n_source+12):(48*n_source+23), 'sub_target'] = follower

            # Follower to leader
            temp_df.loc[(48*n_source+24):(48*n_source+35), 'ste'] = [symb_transfer_entropy(follower_trial_s[n_source, :], leader_trial_t[i, :], delay, 3) for i in range(leader_trial_t.shape[0])]
            temp_df.loc[(48*n_source+24):(48*n_source+35), 'sub_source'] = follower
            temp_df.loc[(48*n_source+24):(48*n_source+35), 'sub_target'] = leader

            # Follower to follower
            temp_df.loc[(48*n_source+36):(48*n_source+47), 'ste'] = [symb_transfer_entropy(follower_trial_s[n_source, :], follower_trial_t[i, :], delay, 3) for i in range(leader_trial_t.shape[0])]
            temp_df.loc[(48*n_source+36):(48*n_source+47), 'sub_source'] = follower
            temp_df.loc[(48*n_source+36):(48*n_source+47), 'sub_target'] = follower

        print("--- %i Done in %s seconds ---" % ((n_trial + 1), (time.time() - start_time)))
        if n_trial == 0:
            grand_df = temp_df
        else:
            grand_df = grand_df.append(temp_df)
        grand_df.to_csv('/home/horozco/2-STE_graham/csv_grand_matrices/' + pair + '_' + freq_band_s + '_' + freq_band_t + '_' + str(delay) + '.csv')

main()
