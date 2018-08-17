""" Normalized symbolic transfer entropy

This module includes all the necessary functions to calculate normalized
symbolic transfer entropy using NumPy and Matplotlib for the project
Hypermusic, specifically the movement data.

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

import glob
import xlrd

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd




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

# This function intakes the xls and outputs the STE score
def xls_to_ste(leaderFile, followerFile, dim, tau, delay):
    # Initialize parameters
    workbook = xlrd.open_workbook(leaderFile)
    worksheet = workbook.sheet_by_index(0)
    source = worksheet.col_values(3)
    source = np.asarray(source[3:])
    workbook = xlrd.open_workbook(followerFile)
    worksheet = workbook.sheet_by_index(0)
    target = worksheet.col_values(3)
    target = np.asarray(target[3:])

    # Correct mistakes
    source = source[~np.isnan(source)]
    target = target[~np.isnan(target)]
    min_length = min(len(source), len(target))
    source = source[:min_length]
    target = target[:min_length]

    # Take four samples and average them together
    temp_source = np.zeros(int(np.floor(source.shape[0] / 4.)))
    y = 0
    for x in range(1, source.shape[0]):
        if x % 4 == 0:
            temp_source[y] = np.mean((source[x], source[x - 1],
                                     source[x-2], source[x-3]))
            y += 1

    temp_target = np.zeros(int(np.floor(target.shape[0] / 4.)))
    y = 0
    for x in range(1, target.shape[0]):
        if x % 4 == 0:
            temp_target[y] = np.mean((target[x], target[x - 1],
                                      target[x - 2], target[x - 3]))
            y += 1

    # Get z score for each the trial
    source_z = sp.stats.mstats.zscore(temp_source)
    target_z = sp.stats.mstats.zscore(temp_target)


    # Calculate ste
    data = np.concatenate((np.expand_dims(target_z, axis=1), #follower, target
                           np.expand_dims(source_z, axis=1)), #leader, source
                          axis=1)
    temp = NSTE(data, dim, tau, delay, False)
    return temp

#######################################################
##                  From XLS to STE                  ##
#######################################################
# A few notes...
# -> For first iteration of analysis, I decided to only use a tau of 1 and
#    a delay of 1!
# ->Remember: the NSTE function outputs a STE list which is interpreted as:
#     ->Column 1: source to target (y to x)
#     ->Column 2: target to source (x to y)

# Please note that I removed manually the dummy files (first recording of
# each thing) and the "doubled" file from P02 (we recorded twice the
# first one).

# Initialize parameters
for name in glob.glob('/Users/hectorOrozco/Desktop/hm_nuri/*.xls'):
    filenames_xls.append(name)
leaderFile_h = [ '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hMelA_video3_frames_0_10049_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hMelA_video5_frames_0_10403_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hValA_video1_frames_0_12803_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hValA_video3_frames_0_11687_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hMelB_video8_frames_0_10049_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hMelB_video10_frames_0_9273_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hValB_video6_frames_0_12083_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hValB_video8_frames_0_12120_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_hMelA_video8_frames_0_4422_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_hMelA_video10_frames_0_4571_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_hValA_video4_frames_0_4911_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_hValA_video6_frames_0_4942_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_hMelB_video3_frames_0_4394_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_hMelB_video5_frames_0_3945_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_hValB_video7_frames_0_5294_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_hValB_video9_frames_0_5135_flow_signals.xls',
    ] # please note that I sorted these manually
followerFile_h = [ '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hMelA_video4_frames_0_10049_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hMelA_video6_frames_0_10360_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hValA_video2_frames_0_12855_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_hValA_video4_frames_0_11687_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hMelB_video7_frames_0_10049_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hMelB_video9_frames_0_9367_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hValB_video5_frames_0_12083_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_hValB_video7_frames_0_12224_flow_signals.xls',

 '/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_hMelA_Video7_frames_0_4524_flow_signals.xls',
                 '/Users/hectorOrozco/Desktop/hm_nuri'
                 '/hM_4B_hMelA_video9_frames_0_4571_flow_signals.xls',
                 '/Users/hectorOrozco/Desktop/hm_nuri'
                 '/hM_4B_hValA_video3_frames_0_4969_flow_signals.xls',
                 '/Users/hectorOrozco/Desktop/hm_nuri'
                 '/hM_4B_hValA_video5_frames_0_4890_flow_signals.xls',
                 '/Users/hectorOrozco/Desktop/hm_nuri'
                 
                 '/hM_3A_hMelB_video4_frames_0_4440_flow_signals.xls',
                 '/Users/hectorOrozco/Desktop/hm_nuri'
                 '/hM_3A_hMelB_video6_frames_0_4050_flow_signals.xls',
                 '/Users/hectorOrozco/Desktop/hm_nuri'
                 '/hM_3A_hValB_video8_frames_0_5508_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_hValB_video10_frames_0_5222_flow_signals.xls'
    ]
leaderFile_p = [
'/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCanA_video5_frames_0_7610_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCanA_video7_frames_0_7968_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCloA_video5_frames_0_9324_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCloA_video7_frames_0_9496_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCanB_video2_frames_0_8158_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCanB_video4_frames_0_7907_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCloB_video2_frames_0_11140_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCloB_video4_frames_0_9729_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCanA_video10_frames_0_4742_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCanA_video8_frames_0_3807_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCloA_video4_frames_0_4156_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCloA_video6_frames_0_4363_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCanB_video3_frames_0_3487_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCanB_video5_frames_0_3388_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCloB_video7_frames_0_4174_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCloB_Video9_frames_0_3885_flow_signals.xls'
    ]
followerFile_p = [
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCanA_video6_frames_0_7724_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCanA_video8_frames_0_7968_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCloA_video6_frames_0_9455_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2B_pCloA_video8_frames_0_9583_flow_signals.xls',

 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCanB_video1_frames_0_8158_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCanB_video3_frames_0_7907_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCloB_video1_frames_0_11140_flow_signals.xls',
 '/Users/hectorOrozco/Desktop/hm_nuri/hM_2A_pCloB_video3_frames_0_10027_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCanA_video9_frames_0_4742_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCanA_video7_frames_0_3807_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCloA_video3_frames_0_4137_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_4B_pCloA_video5_frames_0_4329_flow_signals.xls',

'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCanB_video4_frames_0_3524_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCanB_video6_frames_0_3388_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCloB_video8_frames_0_4023_flow_signals.xls',
'/Users/hectorOrozco/Desktop/hm_nuri/hM_3A_pCloB_video10_frames_0_3946_flow_signals.xls',
    ]
dim = 3 # for these first analysis I use very basic parameters
tau = 1
delay = 1
movement_ste_h = np.zeros((len(leaderFile_h), 2)) # Output matrices
movement_ste_p = np.zeros((len(leaderFile_h), 2))

# Loop across trials and get STE separately for homophonic and polyphonic
for z in range(0, len(leaderFile_h)):
    temp = xls_to_ste(leaderFile_h[z], followerFile_h[z], dim, tau, delay)
    movement_ste_h[z, :] = temp # L2F, F2L

for z in range(0, len(leaderFile_p)):
    temp = xls_to_ste(leaderFile_p[z], followerFile_p[z], dim, tau, delay)
    movement_ste_p[z, :] = temp # L2F, F2L


# Transform data into pandas data frame for more easy plotting


df_1 = pd.DataFrame({'Symbolic Transfer Entropy [bits]': movement_ste_h[:,0],
                     'Direction': ['Leader to Follower'] * 16,
                     'Duo Type':['Homophonic'] * 16})

df_2 = pd.DataFrame({'Symbolic Transfer Entropy [bits]': movement_ste_h[:,1],
                     'Direction': ['Follower to Leader'] * 16,
                     'Duo Type':['Homophonic'] * 16})

df_3 = pd.DataFrame({'Symbolic Transfer Entropy [bits]': movement_ste_p[:,0],
                     'Direction': ['Leader to Follower'] * 16,
                     'Duo Type':['Polyphonic'] * 16})

df_4 = pd.DataFrame({'Symbolic Transfer Entropy [bits]': movement_ste_p[:,1],
                     'Direction': ['Follower to Leader'] * 16,
                     'Duo Type':['Polyphonic'] * 16})


movement_ste = df_1.append([df_2, df_3, df_4])

movement_ste['Duo Type'] = movement_ste['Duo Type'].astype('category')
movement_ste['Direction'] = movement_ste['Direction'].astype('category')



plt.clf()
sns.set_style('darkgrid', {'legend.frameon': True})
g = sns.boxplot(x="Duo Type", y='Symbolic Transfer Entropy [bits]',
            hue="Direction", data=movement_ste, palette="Set3")
g.figure.savefig('/Users/hectorOrozco/Desktop/hm_nuri/movement_ste_horizontal.eps'
                 , format='eps', dpi=1000)



# From Python to R
movement_ste.to_csv('/Users/hectorOrozco/Desktop/hm_nuri/movement_ste_horizontal.csv')