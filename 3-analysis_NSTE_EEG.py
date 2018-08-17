""" Normalized symbolic transfer entropy

This module includes all the necessary functions to calculate normalized
symbolic transfer entropy using NumPy and Matplotlib for the project
Hypermusic.

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
Last updated: 08/03/18
"""

import os

import matplotlib.pyplot as plt
import mne
import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd



####################################################################
##                   Filter Consruction Module                    ##
####################################################################
def kaiser_frequency_bands(sfreq, width, ripple_db, figureOut):
    '''
    This function takes in three parameters and outputs 5 FIR filters (wind.
    sinc function using a Kaiswer window) that will go from 1 to 50 Hz to
    isolate frequency bands in EEG. Kaiswer window: concentrates all power
    in main lobe. The parameters we were using: 2400.00 Frequency rate, 2.5
    Hz width, and 60. Hertz as ripple_dB in the stopband
    Args:
        sfreq (float): Sampling frequency of data
        width (float): Width in Hz of the transition band
        ripple_db (float): Ripple in the stopband
        figureOut (boolean): If True, output figure
    Returns:
        filters: (list): a list containing the frequency band name, the freq
        borders, the length of the filter, and the coefficients
    '''
    iter_freqs = [
            ('Delta', 1., 3., '#E45F56'),
            ('Theta', 4., 7., '#A3D39C'),
            ('Alpha', 8., 12., '#7ACCC8'),
            ('Beta', 13., 29., '#4AAAA5'),
            ('Gamma', 30., 45., '#35404F')]
    filters = []
    nyq_rate = sfreq / 2.
    width = width / nyq_rate  # Hz / nyq_rate
    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for band, fmin, fmax, fcolor in iter_freqs:
        N, beta = sp.signal.kaiserord(ripple_db, width)
        taps = sp.signal.firwin(N, [fmin, fmax],
                                window=('kaiser', beta), nyq = nyq_rate,
                                pass_zero = False)
        filters.append((band, fmin, fmax, N, taps))
        if figureOut:
            w, h = sp.signal.freqz(taps, worN = 16000)
            x = (w / np.pi) * nyq_rate
            y = 20 * np.log10(abs(h) ** 2) #Because we are applying filtfilt
            plt.plot(x[0:670], y[0:670], color = fcolor)
            plt.ylabel('Amplitude (dB)')
            plt.xlabel('Frequency (Hz)')
            angles = np.unwrap(np.angle(h))
            plt.plot(x[0:670], angles[0:670], color = fcolor)
            plt.grid()
    if figureOut:
        plt.title('Digital filter frequency response (filtfilt output)')
        ax2 = ax1.twinx()
        plt.ylabel('Angle (radians)')
        plt.axis('tight')
        plt.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                'figs/6-filter_kaiser.eps')
    return filters


######################################################################
##                       Time Frequency Module                      ##
######################################################################
def baseline_values(lBaselineFile, fBaselineFile, generalPath):
    '''
    This function intakes the path of the .set eeglab files to calculate
    the baseline values for (either) a dB or a % change normalization (see
    Cohen's Analyzing Neural Time series Data: Theory and Practice book,
    chapter 18 for a discussion on this
    Args:
        lBaselineFile (string): name of the baseline file of the leader
        fBaselineFile (string): path where all the above files are found
        generalPath (string): Path to where all files are located
    Returns:
        lBaseline (list): list of floats corresponding to the leader's
        frequency band to use either on a percentage normalization or
        a dB norm
        fBaseline (list): list of floats corresponding to the follower's
        frequency band to use either on a percentage normalization or
        a dB norm
        transition_widths (list): list of the transition widths used for each
        filter
    '''
    # Initialize parameters
    iter_freqs = [ ('Delta', 1., 3.), ('Theta', 4., 7.), ('Alpha', 8., 12.),
                   ('Beta', 13., 29.), ('Gamma', 30., 45.)]
    transition_widths = []
    lBaseline = []
    fBaseline = []
    temp = mne.io.read_raw_eeglab(generalPath + lBaselineFile,
                                   preload=False)
    list_channels = [x.encode('utf-8') for x in temp.ch_names]
    list_channels = list_channels[:-1]

    # Filter, get envelope (HT), and get mean activity per freqband and channel
    # Note: Blackman window used to minimize spectral leakage
    for band, fmin, fmax in iter_freqs:
        transition_widths.append([min(max(fmin * 0.25, 2.), fmin),
                                  min(max(fmax * 0.25, 2.), 2400 / 2. - fmax)])
        l_baseline = mne.io.read_raw_eeglab(generalPath + lBaselineFile,
                                            preload=True)

        l_baseline.pick_types(eeg=True, stim=False)
        l_baseline.filter(fmin, fmax, n_jobs = 1,
                          fir_design = 'firwin', method = 'fir',
                          phase = 'zero-double',
                          fir_window = 'blackman')
        l_baseline.apply_hilbert(envelope = True)
        l_baseline._data = l_baseline._data ** 2 # Calculate power
        l_baseline.plot(show = False, scalings = {'eeg': 20e-12})
        plt.title('Leader: {}'.format(band))
        plt.show()
        lBaseline.append(((band, fmin, fmax),
                          np.mean(l_baseline._data, axis = 1)))

        f_baseline = mne.io.read_raw_eeglab(generalPath + fBaselineFile,
                                            preload=True)
        f_baseline.pick_types(eeg=True, stim=False)
        f_baseline.filter(fmin, fmax, n_jobs=1,
                          fir_design='firwin', method='fir',
                          phase='zero-double',
                          fir_window='blackman')
        f_baseline.apply_hilbert(envelope=True)
        f_baseline._data = f_baseline._data ** 2  # Calculate power
        f_baseline.plot(show=False, scalings={'eeg': 20e-12})
        plt.title('Follower: {}'.format(band))
        plt.show()
        fBaseline.append(((band, fmin, fmax),
                          np.mean(f_baseline._data, axis=1)))

    return lBaseline, fBaseline, list_channels
def EEG_to_normTF(leaderFile, followerFile, generalPath, lBaseline, fBaseline):
    '''
    This function intakes the paths to the .set files of both the leader
    and the follower to process them by filtering the individual frequency
    bands and normalizing using a % change from baseline
    Args:
        leaderFile (string): name of the leader's EEG file to be analyzed
        followerFile (string): name of the follower's EEG file to be analyzed
        generalPath (string): Path to where all files are located
        filts (list): a list with the filtering coefficients to use
        lBaseline (numpy.array): list of floats corresponding to the leader's
        frequency band to use either on a percentage normalization or
        a dB norm
        fBaseline (numpy.array): list of floats corresponding to the follower's
        frequency band to use either on a percentage normalization or
        a dB norm
    Returns:
        l_freq_bands (list): list containing the name of the frequency band,
        the frequency limits, and an MNE RawEEGLAB object with the activity at
        the particular frequency band expressed as % change from  baseline
        l_freq_bands (list): list containing the name of the frequency band,
        the frequency limits, and an MNE RawEEGLAB object with the activity at
        the particular frequency band expressed as % change from  baseline
    '''
    # Initialize parameters
    iter_freqs = [ ('Delta', 1., 3.), ('Theta', 4., 7.), ('Alpha', 8., 12.),
                   ('Beta', 13., 29.), ('Gamma', 30., 45.)]
    l_freq_bands = []
    f_freq_bands = []
    x = 0

    # Filter, get envelope (HT), and get mean activity per freqband and channel
    # and normalize using baseline values.
    # Note: Blackman window used to minimize spectral leakage
    for band, fmin, fmax in iter_freqs:
        # Leader
        l_EEG = mne.io.read_raw_eeglab(generalPath + leaderFile,
                                            preload=True)
        l_EEG.pick_types(eeg=True, stim=False)
        l_EEG.filter(fmin, fmax, n_jobs = 1, fir_design = 'firwin',
                     method = 'fir',phase = 'zero-double',
                     fir_window = 'blackman')
        l_EEG.apply_hilbert(envelope = True)
        l_EEG._data = l_EEG._data ** 2  # Calculate power
        leader_baseline = np.repeat(np.expand_dims(lBaseline[x][1], axis=1),
                                    l_EEG._data.shape[1], axis=1)
        l_EEG._data = 100 * ((l_EEG._data - leader_baseline)
                             / leader_baseline)
        l_EEG.plot(show = False, scalings = {'eeg': 500},)
        plt.title('Leader: {}'.format(band))
        plt.show()
        l_freq_bands.append(((band, fmin, fmax), l_EEG))

        # Follower
        f_EEG = mne.io.read_raw_eeglab(generalPath + followerFile,
                                       preload=True)
        f_EEG.pick_types(eeg=True, stim=False)
        f_EEG.filter(fmin, fmax, n_jobs=1, fir_design='firwin',
                     method='fir', phase='zero-double',
                     fir_window='blackman')
        f_EEG.apply_hilbert(envelope=True)
        f_EEG._data = f_EEG._data ** 2  # Calculate power
        follower_baseline = np.repeat(np.expand_dims(fBaseline[x][1], axis=1),
                                    f_EEG._data.shape[1], axis=1)
        f_EEG._data = 100 * ((f_EEG._data - follower_baseline)
                             / follower_baseline)
        f_EEG.plot(show=False, scalings={'eeg': 500}, )
        plt.title('Follower: {}'.format(band))
        plt.show()
        f_freq_bands.append(((band, fmin, fmax), f_EEG))
        x += 1
    return l_freq_bands, f_freq_bands


####################################################################
##                          NSTE Module                           ##
####################################################################
#Still gotta do the STE statistical significance function!
def prediction_time(x,y,maxlag):
    '''Select the prediction time in STE as the time lag for the maximum
    cross-correlation between two signals. If it is 0,
    then set it to 1, otherwise, use the time lag. Delta is basically
    how much you delay the target and source "past" from the "future", or the
    points you want to predict
    Args:
        x (numpy.ndarray): the target to be predicted
        y (numpy.ndarray): the source of the prediction
        maxlag (int): number of lags to show
    Returns:
        pred_time (int): the lag that yields the max correlation between
        signals (in samples)
    '''
    (lags, c, line, b) = plt.xcorr(x, y, maxlags = maxlag)
    ik = np.argmax(c)
    tlag_max_cc = lags[ik]
    if tlag_max_cc <= 1:
        pred_time = 1
    else:
        pred_time = tlag_max_cc
    return pred_time
def delay_recons(data, lag, dim):
    '''
    Wrangles data to make analysis easier. Takes in a 2D matrix
    (time, signals) and outputs a 3D matrix
    (time, embedded dimensions, lagged signals). Each column will be lagged
    based on tau
    Args:
        data (numpy.ndarray): data in the form of times, [target, source]
        lag (int or float) :tau, the lag between the source and the target
        dim (int): embedding dimension for the symbolization process
    Returns:
        delayed_data (numpy.ndarray): 3D array (time, embedded
        dimensions, signals) with all the embedded dimensions in one place
        (makes it easier calculate probabilities this way)
    '''
    MaxEpoch = data.shape[0]
    signal = data.shape[1]
    delayed_data = np.zeros((MaxEpoch - lag*(dim-1), dim, signal))
    for c in range(0, signal):
        for j in range(1, dim+1):
            delayed_data[:, j - 1, c] = data[(j-1)*lag:
                                             MaxEpoch - (dim-j)*lag, c]
    return delayed_data
def symb_to_integer(delayed_data_symb):
    '''
    Takes the reconstructed delayed matrix and transform it into
    an Integer column matrix that has integer numbers corresponding to the
    symbolization process. With this matrix, probabilities are calculated
    easily. This is achieved by ranking each of the row's columns and then
    concatenating those numbers using 10 and its powers (i.e. if the ranking
    is 1 3 2, then 1 * 10**2 + 3 * 10**1 + 2 * 10**0 = 132)
    Args:
        delayed_data_symb (numpy.ndarray): 2D matrix of data, embedding
        dimension (the output of delay_recons of either the target or the
        source signal)
    Returns:
        Integer (numpy.ndarray): An array of data with the "symbols" (the
         ranking of the data) encoded as integers.
    '''
    symbols = np.argsort(delayed_data_symb, axis = 1)
    SZ = symbols.shape
    Integer = np.zeros((SZ[0], 1))
    for i in range(0, SZ[1]):
        Integer = Integer + np.expand_dims((symbols[:, i] + 1)
                                           * 10**(SZ[1] - i - 1), 1)
    Integer = Integer.astype(int)
    return Integer
def estimate_prob(Integer_prob):
    '''In-module function that takes in an Integer array (e.g. {X(n + delta),
    X(n), Y(n)} or {X(n), Y(n)} and outputs the probability for each row
    (combination of integers).
    Args:
        Integer_prob (numpy.ndarray, int): the Integer array to calculate
        the probabilities on
    Returns:
        prob (numpy.ndarray, probabilities): an array of the same length
        as the input with the probabilities associated with each of the
        row's combination of symbols/integers
    '''
    SZ = [Integer_prob.shape[0], 1]
    count = np.zeros(SZ)
    b = np.unique(Integer_prob) #discrete/non-overlapping events in INT input
    for i in range(0, b.shape[0]):
        index = np.where(b[i] == Integer_prob) #get the indexes of events
        count[index] = len(index[0]) #number of events
    prob = count / float(SZ[0]) #get probability of all events
    return prob
def integer_to_prob(Integer1, Integer2, Integer3):
    '''Take the incoming integer matrices (output from symb_to_integer) and
    calculate the probabilites needed for NSTE.
    Args:
        Integer1 (numpy.ndarray, int): The "integer" time series that
        corresponds to X(n+delta). Basically, the "Target future"
        Integer2 (numpy.ndarray, int): The "integer" time series that
        corresponds to X(n). Basically, the "Target present"
        Integer3 (numpy.ndarray, int): The "integer" time series that
        corresponds to Y(n). Basically, the "Source present"
    Returns:
        P1 (numpy.ndarray, probabilities): P[X(n + delta), X(n), Y(n)]
        P2 (numpy.ndarray, probabilities): P[X(n), Y(n)]
        P3 (numpy.ndarray, probabilities): P[X(n) + delta), X(n)]
        P4 (numpy.ndarray, probabilities): P[X(n)]
    '''

    # Infer the embedding dimension from the data structures
    SymbolLen = max(np.ceil(np.log10(Integer1 + 0.1)))
    SymbolLen = max(SymbolLen, max(np.ceil(np.log10(Integer2 + .1))))
    SymbolLen = max(SymbolLen, max(np.ceil(np.log10(Integer3 + .1))))

    # INT1 = {X(n + delta), X(n), Y(n)}; {Target future, Target past,
    # Source past}
    INT1 = Integer1 * 10**(SymbolLen * 2) \
           + Integer2 * 10**(SymbolLen * 1) \
           + Integer3
    INT1 = INT1.astype(int)

    # INT2 = {X(n), Y(n)}; {Target past, Source past}
    INT2 = Integer2 * 10**(SymbolLen * 1) + Integer3
    INT2 = INT2.astype(int)

    # {X(n + delta), X(n)}; {Target future, Target past}
    INT3 = Integer1 * 10**(SymbolLen * 1) + Integer2
    INT3 = INT3.astype(int)

    # {X(n)}; {Target past}
    INT4 = Integer2

    P1 = estimate_prob(INT1) # P[X(n) + delta), X(n), Y(n)]
    P2 = estimate_prob(INT2) # P[X(n), Y(n)]
    P3 = estimate_prob(INT3) # P[X(n) + delta), X(n)]
    P4 = estimate_prob(INT4) # P[X(n)]

    # P1 gives us the set of possible outcomes
    (temp, U_Index) = np.unique(INT1, return_index = True)
    P1 = P1[U_Index]
    P2 = P2[U_Index]
    P3 = P3[U_Index]
    P4 = P4[U_Index]
    return P1, P2, P3, P4
def NSTE(data, dim, lag, delta, norm):
    '''
    Calculates the normalized symbolic transfer entropy of the two signals
    in data. It only applies to bivariate data. The calculation is done
    column (the whole "set" of possible events) to row (specific event)
    Args:
        data (numpy.ndarray): data in the form of times, [target, source].
        Note that is has to be an [x, 2] array
        dim (int): embedding dimension for the symbolization process
        lag (int or float): tau, the lag between the source and the target
        delta (int): the lag between the signal and the prediction
        norm (boolean): calculate normalization or not
    Returns:
        symb_trans_entropy (numpy.ndarray): a [# of taus, 2] array with the
        raw symbolic transfer entropy values. Column 1: source to target
        (y to x). Column 2: target to source (x to y)
        norm_symb_trans_entropy (numpy.ndarray): a [# of taus, 2] array with
        the normalized symbolic transfer entropy values. Column 1: shuffled
        source to target (y to x). Column 2: target to shuffled source (x to y)
    '''
    ch = data.shape[1]
    ###############################
    # Part1. STE of Original Data #
    ###############################
    Ddata = delay_recons(data, lag, dim) # Wrangle data to make analyses easier
    INT = np.zeros((len(Ddata), ch), dtype = int)
    for c in range(0, ch):
        INT[:, c] = np.squeeze(symb_to_integer(Ddata[:,:,c]))
    Int_future = INT[delta:, :]
    Int_past = INT[0:-delta,:]
    #Actually compute the probability
    (P1, P2, P3, P4) = integer_to_prob(
            Int_future[:, 0], Int_past[:, 0], Int_past[:, 1]
            )
    symb_trans_entropy = np.empty([2])
    norm_symb_trans_entropy = np.empty([2])

    # Note: P(A| B, C) = P(A and B and C)/ P(B and C)
    symb_trans_entropy[0] = sum(P1 * (np.log2(P1 * P4)
                     - np.log2(P2 * P3)))
    (P1, P2, P3, P4) = integer_to_prob(
            Int_future[:, 1], Int_past[:, 1], Int_past[:, 0]
            )
    symb_trans_entropy[1] = sum(P1 * (np.log2(P1 * P4)
                     - np.log2(P2 * P3)))

    if norm:
        ###############################
        # Part2. STE of Shuffled Data #
        ###############################
        # Shuffle Y (source) to calculate the "bias" and normalize STE
        data2 = data[:,1]
        length = data2.shape[0]
        halflen = np.floor(length/2).astype(int)
        shuffled_data2 = np.concatenate((data2[halflen:], data2[0:halflen]), axis=0)
        Ddata = delay_recons(np.stack((data[:,0], shuffled_data2), axis = 1), lag, dim)
        INT = np.zeros((Ddata.shape[0], ch), dtype=int)
        for c in range(0, ch):
            INT[:, c] = np.squeeze(symb_to_integer(Ddata[:,:,c]))
        Int_future = INT[delta:, :]
        Int_past = INT[0:-delta,:]
        (P1, P2, P3, P4) = integer_to_prob(
                Int_future[:, 0], Int_past[:, 0], Int_past[:, 1]
                )
        sft_T = sum(P1 * (np.log2(P1 * P4)
                         - np.log2(P2 * P3)))
        H_XY = sum(P1 * (np.log2(P2)
                         - np.log2(P1 * P2)))
        norm_symb_trans_entropy[0] = (symb_trans_entropy[0] - sft_T)/H_XY

        ###############################
        # Part3. STE of Shuffled Data #
        ###############################
        # Shuffle X (target) to calculate the "bias" and normalize STE
        data1 = data[:, 0]
        length = data1.shape[0]
        shuffled_data1 = np.concatenate((data1[halflen:], data1[0:halflen]),
                                        axis=0)
        Ddata = delay_recons(np.stack((shuffled_data1, data[:, 1]), axis=1), lag,
                             dim)
        INT = np.zeros((Ddata.shape[0], ch), dtype=int)
        for c in range(0, ch):
            INT[:, c] = np.squeeze(symb_to_integer(Ddata[:, :, c]))
        Int_future = INT[delta:, :]
        Int_past = INT[0:-delta, :]

        (P1, P2, P3, P4) = integer_to_prob(
                Int_future[:, 1], Int_past[:, 1], Int_past[:, 0]
                )
        sft_T = sum(P1 * (np.log2(P1 * P4)
                          - np.log2(P2 * P3)))
        H_XY = sum(P1 * (np.log2(P2)
                         - np.log2(P1 * P2)))
        norm_symb_trans_entropy[1] = (symb_trans_entropy[1] - sft_T) / H_XY
        return symb_trans_entropy, norm_symb_trans_entropy
    return symb_trans_entropy
def NSTE_oneWay(data, dim, lag, delta):
    '''
    Calculates the normalized symbolic transfer entropy of the two signals
    in data. It only applies to bivariate data. The calculation is done
    column (the whole "set" of possible events) to row (specific event)
    Args:
        data (numpy.ndarray): data in the form of times, [target, source].
        Note that is has to be an [x, 2] array
        dim (int): embedding dimension for the symbolization process
        lag (int or float): tau, the lag between the source and the target
        delta (int): the lag between the signal and the prediction
        norm (boolean): calculate normalization or not
    Returns:
        symb_trans_entropy (numpy.ndarray): a [# of taus, 2] array with the
        raw symbolic transfer entropy values. Column 1: source to target
        (y to x). Column 2: target to source (x to y)
        norm_symb_trans_entropy (numpy.ndarray): a [# of taus, 2] array with
        the normalized symbolic transfer entropy values. Column 1: shuffled
        source to target (y to x). Column 2: target to shuffled source (x to y)
    '''
    ch = data.shape[1]
    ###############################
    # Part1. STE of Original Data #
    ###############################
    Ddata = delay_recons(data, lag, dim) # Wrangle data to make analyses easier
    INT = np.zeros((len(Ddata), ch), dtype = int)
    for c in range(0, ch):
        INT[:, c] = np.squeeze(symb_to_integer(Ddata[:,:,c]))
    Int_future = INT[delta:, :]
    Int_past = INT[0:-delta,:]
    #Actually compute the probability
    (P1, P2, P3, P4) = integer_to_prob(
            Int_future[:, 0], Int_past[:, 0], Int_past[:, 1]
            )
    symb_trans_entropy = sum(P1 * (np.log2(P1 * P4)
                     - np.log2(P2 * P3)))
    return symb_trans_entropy
def nste_statistical_significance():
    NumWin = 1  #
    winSize = np.around((frequency_map[0][1].data.shape[1]) / NumWin)
    TotalWin = np.floor(
        (frequency_map[0][1].data.shape[1]) / winSize)  #Sanity check
    # Yes, we still get 5 total windows
    RanWin = np.random.permutation(int(TotalWin))
    UsedWin = np.sort(RanWin)
    p_value = UsedWin
    return p_value


####################################################################
##                  From EEG to Channel-wise STE                  ##
####################################################################
# A few notes...
# -> For first iteration of analysis, I decided to only use a tau of 1 and
#    a delay of 1!
# ->Remember: the NSTE function outputs a STE list which is interpreted as:
#     ->Column 1: source to target (y to x)
#     ->Column 2: target to source (x to y)

# Initialize parameters
leaderFile =  'M_4B_pcloB_EEG42017.12.01_17.08.27.set'
followerFile = 'M_3A_pclob_eeg42017.12.01_17.17.06.set'


ch_groupings = ['Frontal', 'Frontal', 'Frontal', 'Frontal', 'Frontal',
                'Frontal', 'Frontal', 'Frontal', 'Frontal', 'Frontal',
                'Frontal', 'Frontal', 'Frontal', 'Frontal', 'Frontal',
                'Frontal',
                'Temporoparietal L', 'Temporoparietal L', 'Temporoparietal L',
                'Central', 'Central', 'Central',
                'Temporoparietal R', 'Temporoparietal R', 'Temporoparietal R',
                'Temporoparietal L', 'Temporoparietal L', 'Temporoparietal L',
                'Central', 'Central', 'Central',
                'Temporoparietal R', 'Temporoparietal R', 'Temporoparietal R',
                'Temporoparietal L', 'Temporoparietal L', 'Temporoparietal L',
                'Central', 'Central', 'Central',
                'Temporoparietal R', 'Temporoparietal R', 'Temporoparietal R',
                'Temporoparietal L', 'Temporoparietal L', 'Temporoparietal L',
                'Parietoccipital', 'Parietoccipital', 'Parietoccipital',
                'Temporoparietal R', 'Temporoparietal R', 'Temporoparietal R',
                'Parietoccipital', 'Parietoccipital', 'Parietoccipital',
                'Parietoccipital', 'Parietoccipital', 'Parietoccipital',
                'Parietoccipital', 'Parietoccipital',
                'Temporoparietal L', 'Temporoparietal R'] #Note that this division is super arbitrary
leaderFile = ['M_3A_hvala_eeg22017.12.01_16.46.09.set',
              'M_3A_hvala_eeg32017.12.01_16.49.51.set',
              'M_3A_hmela_eeg42017.12.01_17.05.29.set',
              'M_3A_hmela_eeg52017.12.01_17.07.33.set',
              'M_3A_pcloa_eeg22017.12.01_17.12.20.set',
              'M_3A_pcloa_eeg32017.12.01_17.13.50.set',
              'M_3A_pcana_eeg42017.12.01_17.26.31.set',
              'M_3A_pcana_eeg52017.12.01_17.28.09.set',
              'M_4B_hvalB_EEG42017.12.01_16.44.01.set',
              'M_4B_hvalB_EEG52017.12.01_16.46.07.set',
              'M_4B_hmelB_EEG22017.12.01_16.52.32.set',
              'M_4B_hmelB_EEG32017.12.01_16.54.33.set',
              'M_4B_pcloB_EEG42017.12.01_17.08.27.set',
              'M_4B_pcloB_EEG52017.12.01_17.10.17.set',
              'M_4B_pcanB_EEG22017.12.01_17.14.24.set',
              'M_4B_pcanB_EEG32017.12.01_17.16.05.set']
followerFile = ['M_4B_hvalA_EEG22017.12.01_16.37.30.set',
                'M_4B_hvalA_EEG32017.12.01_16.41.12.set',
                'M_4B_hmelA_EEG42017.12.01_16.56.52.set',
                'M_4B_hmelA_EEG52017.12.01_16.58.54.set',
                'M_4B_pcloA_EEG22017.12.01_17.03.41.set',
                'M_4B_pcloA_EEG32017.12.01_17.06.00.set',
                'M_4B_pcanA_EEG42017.12.01_17.17.54.set',
                'M_4B_pcanA_EEG52017.12.01_17.19.34.set',
                'M_3A_hvalb_eeg42017.12.01_16.52.40.set',
                'M_3A_hvalb_eeg52017.12.01_16.54.48.set',
                'M_3A_hmelb_eeg22017.12.01_17.01.08.set',
                'M_3A_hmelb_eeg32017.12.01_17.03.15.set',
                'M_3A_pclob_eeg42017.12.01_17.17.06.set',
                'M_3A_pclob_eeg52017.12.01_17.18.48.set',
                'M_3A_pcanb_eeg22017.12.01_17.22.57.set',
                'M_3A_pcanb_eeg32017.12.01_17.24.40.set']

leaderFile = ['M_3A_hvala_eeg22017.12.01_16.46.09.set',
              'M_3A_hmela_eeg42017.12.01_17.05.29.set',
              'M_3A_pcloa_eeg22017.12.01_17.12.20.set',
              'M_3A_pcana_eeg42017.12.01_17.26.31.set',
              'M_4B_hvalB_EEG42017.12.01_16.44.01.set',
              'M_4B_hmelB_EEG22017.12.01_16.52.32.set',
              'M_4B_pcloB_EEG42017.12.01_17.08.27.set',
              'M_4B_pcanB_EEG22017.12.01_17.14.24.set',] #cropped
followerFile = ['M_4B_hvalA_EEG22017.12.01_16.37.30.set',
                'M_4B_hmelA_EEG42017.12.01_16.56.52.set',
                'M_4B_pcloA_EEG22017.12.01_17.03.41.set',
                'M_4B_pcanA_EEG42017.12.01_17.17.54.set',
                'M_3A_hvalb_eeg42017.12.01_16.52.40.set',
                'M_3A_hmelb_eeg22017.12.01_17.01.08.set',
                'M_3A_pclob_eeg42017.12.01_17.17.06.set',
                'M_3A_pcanb_eeg22017.12.01_17.22.57.set',] #cropped

def normTF_to_PDste(generalPath, leaderFile,  l_freq_bands, followerFile, f_freq_bands, dim, tau, delay, list_channels):
    '''

    '''
    # Initialize parameters
    source = l_freq_bands  #leader
    target = f_freq_bands  #follower
    freqs = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']
    total_electrodes = l_freq_bands[0][1]._data.shape[0] +\
                       f_freq_bands[0][1]._data.shape[0]
    ste = np.zeros((total_electrodes, total_electrodes, 5))  #electrodes, electrodes, freq_band

    # Loop through frequency bands
    for w in range(2, 3):
        x = target[w][1]._data #Follower
        y = source[w][1]._data #Leader
        # Leader to follower
        for l_elec_source in range(0, y.shape[0]):
            for f_elec_target in range(0, x.shape[0]):
                crop_int = min(x[f_elec_target].shape, y[l_elec_source].shape)
                data = np.concatenate((np.expand_dims(x[f_elec_target][:crop_int[0]], axis=1),
                                       #follower, x, target
                                       np.expand_dims(y[l_elec_source][:crop_int[0]], axis=1)),
                                      #leader, y, source
                                      axis=1)
                temp = NSTE(data, dim, tau, delay, False)
                ste[l_elec_source, (f_elec_target + 62), w] = temp[0]
                ste[(f_elec_target + 62), l_elec_source, w] = temp[1]

        # Leader to leader
        for l_elec_source in range(0, y.shape[0]):
            for l_elec_target in range(0, y.shape[0]):
                data = np.concatenate((np.expand_dims(y[l_elec_target], axis=1),
                                      #follower, x, target
                                      np.expand_dims(y[l_elec_source], axis=1)),
                                      #leader, y, source
                                      axis=1)
                ste[l_elec_source, l_elec_target, w] = NSTE_oneWay(data, dim, tau, delay)

        # Follower to follower
        for f_elec_source in range(0, x.shape[0]):
            for f_elec_target in range(0, x.shape[0]):
                data = np.concatenate((np.expand_dims(x[f_elec_target], axis=1),
                                      #follower, x, target
                                      np.expand_dims(x[f_elec_source], axis=1)),
                                      #leader, y, source
                                      axis=1)
                ste[(f_elec_source + 62), (f_elec_target + 62), w] = NSTE_oneWay(data, dim, tau, delay)
        ste_pd = pd.DataFrame(ste[:, :, w],
                              index=list_channels + list_channels,
                              columns=list_channels + list_channels)

        # transform STE to a PD array
        ste_pd.to_csv(generalPath + leaderFile[2:4] + '_' + followerFile[2:4] + leaderFile[4:-4] + '_' + freqs[w] + '.csv')


dim = 3  #embedding dimension
tau = 1
delay = 1
generalPath = '/Users/hectorOrozco/Desktop/hM_analysis/output/'

#3A is leading
for file in range(0, (len(leaderFile)/2)):
    lBaselineFile = 'M_3A_baseline2017.12.01_16.35.40.set'
    fBaselineFile = 'M_4B_Baseline2017.12.01_16.27.03.set'
    lBaseline, fBaseline, list_channels = baseline_values(lBaselineFile, fBaselineFile, generalPath)
    l_freq_bands, f_freq_bands = EEG_to_normTF(leaderFile[file], followerFile[file], generalPath, lBaseline, fBaseline)
    normTF_to_PDste(generalPath, leaderFile[file], l_freq_bands, followerFile[file],
                    f_freq_bands, dim, tau, delay, list_channels)

#4B is leading
for file in range(4, len(leaderFile)):
    lBaselineFile = 'M_4B_Baseline2017.12.01_16.27.03.set'
    fBaselineFile = 'M_3A_baseline2017.12.01_16.35.40.set'
    lBaseline, fBaseline, list_channels = baseline_values(lBaselineFile, fBaselineFile, generalPath)
    l_freq_bands, f_freq_bands = EEG_to_normTF(leaderFile[file], followerFile[file], generalPath, lBaseline, fBaseline)
    normTF_to_PDste(generalPath, leaderFile[file], l_freq_bands, followerFile[file],
                    f_freq_bands, dim, tau, delay, list_channels)




# Please note that I manually pasted an copied this!

files_P03 = []

alpha_matrices_files = []
for file in os.listdir(generalPath):
    if file.endswith(".csv"):
        alpha_matrices_files.append(os.path.join(file))



polyphonic_files = ['3A_4B_pcana_eeg42017.12.01_17.26.31_Alpha.csv',
                    '3A_4B_pcloa_eeg22017.12.01_17.12.20_Alpha.csv',
                    '4B_3A_pcanB_EEG22017.12.01_17.14.24_Alpha.csv',
                    '4B_3A_pcloB_EEG42017.12.01_17.08.27_Alpha.csv']

homophonic_files = ['3A_4B_hmela_eeg42017.12.01_17.05.29_Alpha.csv',
                    '3A_4B_hvala_eeg22017.12.01_16.46.09_Alpha.csv',
                    '4B_3A_hmelB_EEG22017.12.01_16.52.32_Alpha.csv',
                    '4B_3A_hvalB_EEG42017.12.01_16.44.01_Alpha.csv',]


# Polyphonic files
df1 = pd.read_csv(generalPath + polyphonic_files[0], header = 0, index_col = 0)
df2 = pd.read_csv(generalPath + polyphonic_files[1], header = 0, index_col = 0)
df3 = pd.read_csv(generalPath + polyphonic_files[2], header = 0, index_col = 0)
df4 = pd.read_csv(generalPath + polyphonic_files[3], header = 0, index_col = 0)
supreme_alpha = df1
for x in range(0, 124):
    for y in range(0,124):
        supreme_alpha.iloc[x][y] = np.mean((df1.iloc[x][y], df2.iloc[x][y],
                                            df3.iloc[x][y], df4.iloc[x][y]))
plt.clf()
sns.set(font_scale = 0.5)
ax = sns.heatmap(supreme_alpha, cmap="Reds", vmax = 0.0002)
ax.invert_yaxis()
plt.yticks(rotation=0)
plt.xticks(rotation='vertical')
ax.figure.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                  'figs/9-polyphonic_alpha.eps', format='eps', dpi=1000)


# Homophonic files
df1 = pd.read_csv(generalPath + homophonic_files[0], header = 0, index_col = 0)
df2 = pd.read_csv(generalPath + homophonic_files[1], header = 0, index_col = 0)
df3 = pd.read_csv(generalPath + homophonic_files[2], header = 0, index_col = 0)
df4 = pd.read_csv(generalPath + homophonic_files[3], header = 0, index_col = 0)
supreme_alpha = df1
for x in range(0, 124):
    for y in range(0,124):
        supreme_alpha.iloc[x][y] = np.mean((df1.iloc[x][y], df2.iloc[x][y],
                                            df3.iloc[x][y], df4.iloc[x][y]))
plt.clf()
sns.set(font_scale = 0.5)
ax = sns.heatmap(supreme_alpha, cmap="Reds", vmax = 0.0002)
ax.invert_yaxis()
plt.yticks(rotation=0)
plt.xticks(rotation='vertical')
ax.figure.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                  'figs/10-homophonic_alpha.eps', format='eps', dpi=1000)




# Correlations between FPz from leader to follower and PMPQ
df1 = pd.read_csv(generalPath + polyphonic_files[0], header = 0, index_col = 0)
df2 = pd.read_csv(generalPath + polyphonic_files[1], header = 0, index_col = 0)
df3 = pd.read_csv(generalPath + polyphonic_files[2], header = 0, index_col = 0)
df4 = pd.read_csv(generalPath + polyphonic_files[3], header = 0, index_col = 0)
df5 = pd.read_csv(generalPath + homophonic_files[0], header = 0, index_col = 0)
df6 = pd.read_csv(generalPath + homophonic_files[1], header = 0, index_col = 0)
df7 = pd.read_csv(generalPath + homophonic_files[2], header = 0, index_col = 0)
df8 = pd.read_csv(generalPath + homophonic_files[3], header = 0, index_col = 0)

x = np.array([df1.iloc[63][1], df2.iloc[63][1], df3.iloc[63][1], df4.iloc[63][1],
             df5.iloc[63][1], df6.iloc[63][1], df7.iloc[63][1], df8.iloc[63][1]
             ])

synergy = ([8, 6, 5.5, 6, 6.5, 7.5, 5.5, 6.5])
synchrony = ([7.5, 4.5, 5, 5.5, 6, 7, 5.5, 6])
quality = ([6.5, 4.5, 5, 5, 6, 7, 4.5, 5])

plt.clf()
data = np.concatenate((np.expand_dims(x, axis=1),
                       np.expand_dims(synergy, axis=1)),
                       axis=1)
data = pd.DataFrame(data = data,
                    columns = ['STE from leader to follower (FCz)',
                               'Synergy (self reported)'])
g = sns.regplot(x ='STE from leader to follower (FCz)',
               y ='Synergy (self reported)', data = data, ci = None,
                color="#FF0000",)
g.set(xlim = (0.00001, 0.00007), ylim = (4, 9))
g.figure.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                  'figs/11-ste_synergy_corr.eps', format='eps', dpi=1000)


plt.clf()
data = np.concatenate((np.expand_dims(x, axis=1),
                       np.expand_dims(synchrony, axis=1)),
                       axis=1)
data = pd.DataFrame(data = data,
                    columns = ['STE from leader to follower (FCz)',
                               'Synchrony (self reported)'])
g = sns.regplot(x ='STE from leader to follower (FCz)',
               y ='Synchrony (self reported)', data = data, ci = None,
                color="#FF0000",)
g.set(xlim = (0.00001, 0.00007), ylim = (4, 9))
g.figure.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                  'figs/12-ste_synchrony_corr.eps', format='eps', dpi=1000)




plt.clf()
data = np.concatenate((np.expand_dims(x, axis=1),
                       np.expand_dims(quality, axis=1)),
                       axis=1)
data = pd.DataFrame(data = data,
                    columns = ['STE from leader to follower (FCz)',
                               'Quality (self reported)'])
g = sns.regplot(x ='STE from leader to follower (FCz)',
               y ='Quality (self reported)', data = data, ci = None,
                color="#FF0000",)
g.set(xlim = (0.00001, 0.00007), ylim = (4, 9))
g.figure.savefig('/Users/hectorOrozco/Desktop/hM_analysis/'
                  'figs/13-ste_quality_corr.eps', format='eps', dpi=1000)
