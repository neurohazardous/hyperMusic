%% Initialize EEGLAB, necessary information...
addpath(genpath('./analysis/_dependencies/eeglab13_4_4b'));
addpath(genpath('./analysis/_dependencies'));
generalpath = '/Users/hectorOrozco/Desktop/hM_analysis/data/P03/eeg';
output = '/Users/hectorOrozco/Desktop/hM_analysis/output/';

%% Iterate through P03 - 3A files
listing = dir('/Users/hectorOrozco/Desktop/hM_analysis/data/P03/eeg/*hM_3A*')
for x = 10:length(listing)
    file = listing(x).name
    chan_locs = '/Users/hectorOrozco/Desktop/hM_analysis/data/P03/eeg/hM_EEGDIGI_3A.sfp';
    EEG = hM_preProEEGLAB(file, generalpath, chan_locs, output);
end

%% Iterate through P03 - 4B files
listing = dir('/Users/hectorOrozco/Desktop/hM_analysis/data/P03/eeg/*hM_4B*')
for x = 1:length(listing)
    file = listing(x).name
    chan_locs = '/Users/hectorOrozco/Desktop/hM_analysis/data/P03/eeg/hM_EEGDIGI_4B.sfp';
    EEG = hM_preProEEGLAB(file, generalpath, chan_locs, output);
end


% Plot the EEG as a sanity check

listing = dir('/Users/hectorOrozco/Desktop/hM_analysis/output/*.set*')
for x = 1:length(listing)
    file = listing(x).name
    EEG = pop_loadset('filename',file,'filepath','/Users/hectorOrozco/Desktop/hM_analysis/output/');
    pop_eegplot( EEG, 1, 1, 1)
    input('Press enter to continue')
    clear EEG
    close all
end
