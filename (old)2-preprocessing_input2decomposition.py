""" Hypermusic Preprocessing

This script intakes EEG data from a gTec system and outputs the preprocessed,
decomposed EEG ready for analysis

The preprocessing steps used are heavily influenced by...
(1) Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., & Robbins, K. A.
(2015). The PREP pipeline: standardized preprocessing for large-scale EEG
analysis. Frontiers in neuroinformatics, 9, 16.

(2) https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline

Written by Hector D Orozco Perez
Last updated: 08/03/18
"""



# HDF stands for "Hierarchical Data Format". Every object in an an HDF5 file has a name, and they are
# arranged in a POSIX-style hierarchy. The folders in the system are called groups. You get access to those with the
# the method h5py.File()
# To access subgroups of data, you can use Python dictionary-style interface using the item-retrieval syntax.
# printname and the visit method to iterate through all the groups folders


# The hdf5 file has a weird XML-like structure were some of the needed data is stored. These XmlDictConfig class has a
# very nice method that is based on the lxml module and takes in a xml and transforms it into a dictionary. Not the most optimal,
# but it allows me to extract the sampling rate quite easily

for name in gTec['AsynchronData/AsynchronSignalTypes']:
    print name

# Needed modules
import mne
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
import matplotlib.pyplot as plt
import scipy as sp



# io_hdf5 takes in a hdf5 gTec file and outputs an mne RawArray. Just specify the file (x, y) = io_hdf5(eegFile)
# The function also needs mne, numpy, h5py, and the xml.etree.ElementTree modules
class XmlListConfig(list):
    def __init__(self, aList):
        for element in aList:
            if element:
                # treat like dict
                if len(element) == 1 or element[0].tag != element[1].tag:
                    self.append(XmlDictConfig(element))
                # treat like list
                elif element[0].tag == element[1].tag:
                    self.append(XmlListConfig(element))
            elif element.text:
                text = element.text.strip()
                if text:
                    self.append(text)


class XmlDictConfig(dict):
    def __init__(self, parent_element):
        if parent_element.items():
            self.update(dict(parent_element.items()))
        for element in parent_element:
            if element:
                # treat like dict - we assume that if the first two tags
                # in a series are different, then they are all different.
                if len(element) == 1 or element[0].tag != element[1].tag:
                    aDict = XmlDictConfig(element)
                # treat like list - we assume that if the first two tags
                # in a series are the same, then the rest are the same.
                else:
                    # here, we put the list in dictionary; the key is the
                    # tag name the list elements all share in common, and
                    # the value is the list itself
                    aDict = {element[0].tag: XmlListConfig(element)}
                # if the tag has attributes, add those to the dict
                if element.items():
                    aDict.update(dict(element.items()))
                self.update({element.tag: aDict})
            # this assumes that if you've got an attribute in a tag,
            # you won't be having any text. This may or may not be a
            # good idea -- time will tell. It works for the way we are
            # currently doing XML configuration files...
            elif element.items():
                self.update({element.tag: dict(element.items())})
            # finally, if there are no child tags and no attributes, extract
            # the text
            else:
                self.update({element.tag: element.text})

def io_hdf5(eegFile):



    # First we open the file (read and write)
    gTec = h5py.File(eegFile, 'r+')

    # There are some XML-like structures that need to be parsed to access certain information (sampling rate...)
    xml = [name.decode() for name in gTec['RawData/AcquisitionTaskDescription']]
    xml = ', '.join(xml)
    root = ElementTree.fromstring(xml)
    xlmdict = XmlDictConfig(root)

    # We get the raw data into a Numpy Array and reshape it (electrodes, samples)
    rawData = np.transpose(np.array(gTec['RawData/Samples']))

    # MNE needs some basic information from this file to create the info object that is then passed to the create raw function
    n_channels = int(xlmdict.get('NumberOfAcquiredChannels'))
    sampling_rate = float(xlmdict.get('SamplingFrequency'))
    ch_names = ['FP1', 'FPz', 'FP2', 'AF7', 'AF3', 'AF4', 'AF8', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2', 'TP9', 'TP10', 'leftEar', 'rightEar']
    ch_types = []
    for x in range(0, len(ch_names)):
        ch_types.append('eeg')
    number_events = len(gTec['AsynchronData/TypeID'])

    #Pass all this info to the create_info function and create the raw object
    info = mne.create_info(ch_names = ch_names, sfreq = sampling_rate, ch_types = ch_types)
    rawEEG = mne.io.RawArray(rawData, info)
    return rawEEG

# function performs a basic preprocessing at this point of things using mne methods
def eeg_preproc(rawEEG):
    # Take out the "problematic channels"
    epochedEEG = rawEEG.pick_types(eeg = True, exclude = ['leftEar', 'rightEar'])
    #Band-pass 1 - 50 Hz
    filteredEEG = epochedEEG.filter(1., 50., h_trans_bandwidth='auto', filter_length='auto', method = 'fir', phase='zero')
    #Import channel general locations
    montage = mne.channels.read_montage('standard_1020')
    chanLocEEG = filteredEEG.set_montage(montage, set_dig = True)
    #Re-reference to common average
    contEEG = chanLocEEG.set_eeg_reference(ref_channels='average', projection=False)

    return contEEG

# function takes in preprocessed EEG MNE data. You need the numpy module, matplotlib.pyplot, and from mne you  need: from mne.time_frequency import tfr_morlet, psd_multitaper
def timeFrequency_bands(contEEG, eegFile):
    gTec = h5py.File(eegFile, 'r+')
    events = np.zeros((1, 3), dtype=np.int)  # Event matrix in MNE always have 3 columns, usually we only care about the first (sample time stamp) and last one (ID of trigger)
    events[0, 0] = int(gTec['AsynchronData/Time'][-1])  # I'm only taking the timing of the last event out of the three square waves that I send
    events[0, 2] = int(gTec['AsynchronData/TypeID'][-1])  # I'm only taking the last event's ID out of the three square waves that I send (they're all the same)
    event_id = dict(music=events[0, 2])
    sampling_rate = float(xlmdict.get('SamplingFrequency'))
    ch_names = ['FP1', 'FPz', 'FP2', 'AF7', 'AF3', 'AF4', 'AF8', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8',
                'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4',
                'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz',
                'P2', 'P4', 'P6', 'P8', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2', 'TP9', 'TP10', 'leftEar',
                'rightEar']
    tmin = 0.0  # because I don't really care about baseline for this experiment, I'm starting at 0 of the last square wave
    tmax = ((contEEG.last_samp - events[0][0])/sampling_rate) - 1 # -1 just to be sure that I am still within the boundaries of the time
    iter_freqs = [
        ('Delta', 1., 3.),
        ('Theta', 4., 7.),
        ('Alpha', 8., 12.),
        ('Beta', 13., 29.),
        ('Gamma', 30., 45.)
    ]
    frequency_map = list()
    for band, fmin, fmax in iter_freqs:
        rawEEG = io_hdf5(eegFile)
        contEEG = eeg_preproc(rawEEG)
        contEEG.filter(fmin, fmax, n_jobs = 1, l_trans_bandwidth = 1, h_trans_bandwidth = 1, fir_design = 'firwin')
        contEEG.apply_hilbert(n_jobs = 1, envelope = False)
        channels = mne.pick_channels(ch_names, include = ['FCz']) #Access FCz
        epochs = mne.Epochs(contEEG, events, event_id, tmin, tmax, picks = channels, baseline = None, reject = None, preload = True)
        epochedEEG = mne.EpochsArray(data=np.abs(epochs.get_data()), info=epochs.info, tmin=epochs.tmin) #You get the absolute value of the complex signal to get the envelope here and epoch the data
        frequency_map.append(((band, fmin, fmax), epochedEEG.average())) #I'm not really averaging anything here because I only have one trial lelz
    return frequency_map #to acces the evokedarray itself, use frequency_map[n][1], where n specifies what frequency band you are interested in

#sampling_rate = float(xlmdict.get('SamplingFrequency'))
#frequency_map[0][1].pick_channels(ch_names = ['FCz']) #get raw data
#    raw = frequency_map[0][1].data









#Plotting the power
fig, axes = plt.subplots(5, 1, figsize=(10, 7), sharex=True, sharey=True)
colors = plt.get_cmap('winter_r')(np.linspace(0, 1, 5))
for ((freq_name, fmin, fmax), average), color, ax in zip(frequency_map, colors, axes.ravel()[::-1]):
    times = average.times * 1e3
    gfp = np.sum(average.data ** 2, axis=0)
    ax.plot(times, gfp, label=freq_name, color=color, linewidth=3)
    ax.axhline(0, linestyle='--', color='grey', linewidth=2)
    ax.grid(True)
    ax.set_ylabel('GFP')
    ax.annotate('%s (%d-%dHz)' % (freq_name, fmin, fmax),
                xy=(0.95, 0.8),
                horizontalalignment='right',
                xycoords='axes fraction')
    ax.set_xlim(tmin, tmax*1000)
    ax.set_ylim(-100, 20000)
axes.ravel()[-1].set_xlabel('Time [ms]')
plt.show()