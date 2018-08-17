function [EEGout] = hM_preProEEGLAB(file, generalpath, chan_locs, output)
%{ 
Hypermusic Preprocessing

Remember: Computations in EEGLAB should be done in double precision, 
specially when running ASR!

This script intakes EEG data from a gTec system and outputs the 
preprocessed, decomposed EEG ready for analysis.

I use EEGLAB and plugins (https://sccn.ucsd.edu/eeglab/index.php) to 
do all the preprocessing. I created this pipeline mainly based on
Makoto's preprocessing pipeline 
(https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline)

Steps: 
1. Load file and delete extra events
2. High pass filter 1 Hz (Hamming windowed sinc FIR)
3. Load channel locations, delete extra channels, and assign 10 20 name
4. Take out line noise using spectral regression
5. Run artifact subspace reconstruction
6. Get info regarding interpolated channels
7. Interpolate bad channels using spherical method
8. Perform CAR 
9. Trim data around trigger event
10. Save output .set file
11. Write bad channels to txt file 

Written by Hector D Orozco Perez
Last updated: 10/03/18
%}

%% Initialize variables
load('./analysis/_dependencies/64ch_withoutA1A2.mat');
ch_names = squeeze(Mon.electrodename);
filename = [file generalpath]

%% Preprocessing
% Load data and delete the first two triggers
EEG = pop_loadhdf5('filename', file, 'filepath', generalpath,... 
                   'rejectchans', [], 'ref_ch', [])
for x = 1:(length(EEG.event)-1)
EEG = pop_editeventvals(EEG,'delete',1);
end               

% High pass filter @ 1Hz (default filtering strategy) using Hamming window
EEG = pop_eegfiltnew(EEG, [], 1, 7920, true, [], 1);
 
% Load channel location file, rename them, and remove earlobes
EEG=pop_chanedit(EEG, 'load',{chan_locs 'filetype' 'sfp'},'changefield',{4 'datachan' 0});
EEG = pop_select( EEG,'nochannel',{'A-63' 'A-64'});
for n = 1:length(ch_names)
   EEG.chanlocs(n).labels = ch_names{n};
end


% Save original EEG before channel rejection/artifact correction
originalChanLocs = EEG.chanlocs;

% Cleanline; note that I used step size == window size (Makoto suggestion)
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:62], ...
                    'computepower',1, 'linefreqs',[60 120], ...
                    'normSpectrum',0,'p',0.01, ...
                    'pad',2,'plotfigures',0,'scanforlines',1, ...
                    'sigtype','Channels','tau',100,'verb',1, ...
                    'winsize',4,'winstep',4);

% Run ASR, use 8 SD because data is very artifactual
% -> high-pass filtering disabled (we are already doing that)
% -> noise based rejection disabled
% -> bad window disabled because we cannot afford to loose data
EEG = clean_rawdata(EEG, 5, [-1], 0.8, -1, 8, -1);
close all

% Get info regarding rejected channels
badChannels = cell(1,length(EEG.etc.clean_channel_mask));
for x = 1:length(EEG.etc.clean_channel_mask)
    if ~EEG.etc.clean_channel_mask(x)
        badChannels{x} = ch_names{x}
    end 
end
badChannels = badChannels(~cellfun('isempty',badChannels)) 

% Interpolate bad channels before doing CAR
EEG = pop_interp(EEG, originalChanLocs, 'spherical');

% Re-reference the data to average by keeping the matrix "full rank"
EEG = fullRankAveRef(EEG);

% Trim data from third event to end of data
triggerTime = EEG.event.latency/EEG.srate
EEGtime = length(EEG.data)/EEG.srate
EEGout = pop_rmdat( EEG, {'Trigger 1'},[-1 (EEGtime-triggerTime)] ,0);

% Save output EEG
EEGout = pop_saveset( EEGout, 'filename',[file(2:end-5) '.set'],'filepath',output);

% Write text file with session bad channels
fid = fopen( './output/badChannels.txt', 'a' );
fprintf(fid, 'File: %s\n', file(2:end-5));
for x=1:length(badChannels)
    if x == 1
        fprintf(fid, 'Interpolated channels: %s, ', badChannels{x});
    elseif x < length(badChannels)  
        fprintf(fid, '%s, ', badChannels{x});
    else 
        fprintf(fid, '%s\n\n', badChannels{x});
    end 
end 
fclose(fid);
clear badChannels

end