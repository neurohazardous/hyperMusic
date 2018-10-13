%% == License ==========================================================
%{
This file is part of the project hyperMusic. All of
hyperMusic code is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
megFingerprinting is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with megFingerprinting.
If not, see <https://www.gnu.org/licenses/>.
%}

%% == Remarks ==========================================================
%{
Remember: Computations in EEGLAB should be done in double precision,
specially when running ASR!

This cortical patch implementation (Limpiti et al., 2006) is based on code
done by Dr Phil Chrapka. You can find the original code here:
https://github.com/pchrapka/fieldtrip-beamforming

This script intakes EEG data from a gTec system and outputs the
preprocessed, decomposed EEG (19 regions of interest) ready for further
analysis.

I use EEGLAB, FieldTrip, and plugins
(https://sccn.ucsd.edu/eeglab/index.php) to
do all the preprocessing. I created this pipeline mainly based on
Makoto's preprocessing pipeline
(https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline)
%}

%% == hyperMusic: Preprocessing =========================================
%{
2) Load file and delete extra events
3) High pass filter 1 Hz (Hamming windowed sinc FIR)
4) Load channel locations and delete extra channels
5) Take out line noise using spectral regression
5) Run artifact subspace reconstruction
6) Get info regarding interpolated channels
7) Interpolate bad channels using spherical method
8) Perform CAR
9) Trim data around trigger event
10) Save output .set file
11) Write bad channels to txt file
1) Headmodel based on based on ICBM 152 atlas
%}

%% == 0) Setup path ======================================================
addpath('dependencies/');
addpath('dependencies/mni152')
addpath('dependencies/spm12/')
addpath('dependencies/eeglab14_1_2b/');
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab
close all
addpath('dependencies/fieldtrip-21092018/')
ft_defaults
global ft_default
ft_default.spmversion='spm12';

%% == 1) Headmodel based on based on ICBM 152 atlas ==================
% This bit is the MRI step of the pipeline
% load MRI data for ICMB 152
mri = ft_read_mri('dependencies/mni152/MNI152_T1_0.5mm.nii');
mri.coordsys = 'mni';

% MNI coordinates are taken from: Cutini S, Scatturin P, Zorzi M (2011): A
% new method based on ICBM152 head surface for probe placement in
% multichannel fNIRS

% check fiducial locations
cfg = [];
cfg.location = [0 84 -43]; % Nasion (NAS)
ft_sourceplot(cfg,mri);

cfg.location = [-75.09 -19.49 -47.98]; % Right pre-auricular point (RPA)
ft_sourceplot(cfg,mri);

cfg.location = [76 -19.45 -47.7]; % Left pre-auricular point (LPA)
ft_sourceplot(cfg,mri);

% add fiducials to mri header
fiducials = [];
fiducials.nas = [0 84 -43];
fiducials.lpa = [-75.09 -19.49 -47.98];
fiducials.rpa = [76 -19.45 -47.7];
fid_names = fieldnames(fiducials);

% transform from MNI to original anatomical coordinate space
transform = mri.transform;
R_inv = inv(transform(1:3,1:3));
d = transform(1:3,4);
inv_transform = [R_inv -R_inv*d; 0 0 0 1];

fid_anatomical = [];
for n_patch=1:length(fid_names)
    field = fid_names{n_patch};
    fid_anatomical.(field) = ft_warp_apply(inv_transform, fiducials.(field), 'homogenous');
    
    % Sanity checks: they should be equal
    fid_mni = ft_warp_apply(transform, fid_anatomical.(field), 'homogenous');
    fprintf('Fiducial: %s\nMNI\n',field);
    disp(fid_mni);
    fprintf('Original\n');
    disp(fiducials.(field));
end

mri.hdr.fiducial.mri = fid_anatomical;
mri.hdr.fiducial.head = fid_anatomical;

% save mri
outputfile = 'dependencies/mni152/icbm152_mri.mat';
save(outputfile,'mri');


% segment volume
% Source: http://www.agricolab.de/template-headmodel-for-fieldtrip-eeg-source-reconstruction-based-on-icbm152/
inputfile = outputfile;
mri = load(inputfile, 'mri');
mri = mri.mri;

cfg = [];
cfg.brainthreshold = 0.5;
cfg.scalpthreshold = 0.15;
cfg.downsample = 1;
cfg.spmversion='spm12';
cfg.output = {'brain','skull','scalp'};

seg = ft_volumesegment(cfg, mri);

outputfile = 'dependencies/mni152/icbm152_seg.mat';
save(outputfile,'seg');


% prepare mesh
inputfile = outputfile;
seg = load(inputfile);
seg = seg.seg;

cfg = [];
cfg.spmversion='spm12';
cfg.method = 'projectmesh';
cfg.tissue = {'brain','skull','scalp'};
cfg.numvertices = [1000, 1000, 1000];

mesh = ft_prepare_mesh(cfg,seg);
scaling = [0.999 1 1.001];
for n_patch=1:length(scaling)
    mesh(n_patch).pos = mesh(n_patch).pos.*scaling(n_patch);
end

outputfile = 'dependencies/mni152/icbm152_mesh.mat';
save(outputfile,'mesh');

figure
hold on
ft_plot_mesh(mesh(1));
ft_plot_mesh(mesh(2));
ft_plot_mesh(mesh(3));
alpha 0.1

% prepare headmodel
inputfile = outputfile;
mesh = load(inputfile);
mesh = mesh.mesh;

cfg = [];
cfg.method = 'dipoli';
cfg.spmversion='spm12';
headmodel = ft_prepare_headmodel(cfg,mesh);
headmodel = ft_convert_units(headmodel,'mm');

outputfile = 'dependencies/mni152/icbm152_bem.mat';
save(outputfile,'headmodel');

% get fiducials
fields = {'nas','lpa','rpa'};
nfields = length(fields);

% allocate mem
fid_pos = zeros(nfields,3);
fid_names = cell(nfields,1);

% get fiducials
transm = mri.transform;
for j=1:nfields
    fid_coord = mri.hdr.fiducial.mri.(fields{j});
    fid_coord = ft_warp_apply(transm, fid_coord, 'homogenous');
    fid_pos(j,:) = fid_coord;
    fid_names{j} = upper(fields{j});
end

% Plot the headmodel just in case!
figure,
ft_plot_mesh(headmodel.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(headmodel.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);
% plot fiducials,
hold on
style = 'or';

% plot points
plot3(fid_pos(:,1), fid_pos(:,2), fid_pos(:,3), style);
% plot labels
for j=1:size(fid_pos,1)
    text(fid_pos(j,1), fid_pos(j,2), fid_pos(j,3), fid_names{j});
end

for n_sub=5:6 %length(subjects)
    %% == 9) Electrode allignment with template ==============================
    % It is easier if you run all the participants at once!
    % Properties
    chan_locs = ['data/raw/' pairs{n_sub} '/eeg/hM_EEGDIGI_' subjects{n_sub} '.sfp'];
    elec = ft_read_sens(chan_locs, 'senstype', 'eeg');
    elec = ft_convert_units(elec, 'mm');
    
    nas=mri.hdr.fiducial.mri.nas;
    lpa=mri.hdr.fiducial.mri.lpa;
    rpa=mri.hdr.fiducial.mri.rpa;
    
    transm=mri.transform;
    
    nas=ft_warp_apply(transm,nas, 'homogenous');
    lpa=ft_warp_apply(transm,lpa, 'homogenous');
    rpa=ft_warp_apply(transm,rpa, 'homogenous');
    
    % plot pre allignment
    figure;
    % head surface (scalp)
    ft_plot_mesh(headmodel.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
    hold on;
    % electrodes
    ft_plot_sens(elec, 'marker', 'o', 'color', 'black');
    
    % create a structure similar to a template set of electrodes to allign
    % MRI with fiducials (first pass)
    fid = [];
    fid.elecpos       = [nas; lpa; rpa];       % ctf-coordinates of fiducials
    fid.label         = {'NZ','LPA','RPA'};    % same labels as in elec
    fid.unit          = 'mm';                  % same units as mri
    
    % First alignment: fiducials
    cfg               = [];
    cfg.method        = 'fiducial';
    cfg.target        = fid;                   % see above
    cfg.elec          = elec;
    cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec
    elec_aligned      = ft_electroderealign(cfg);
    
    % plot post allignment
    figure;
    ft_plot_sens(elec_aligned,'marker', 's', 'color', 'black');
    hold on;
    ft_plot_mesh(headmodel.bnd(1),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp
    
    % we do a second interactive pass just in case!
    cfg           = [];
    cfg.method    = 'interactive';
    cfg.elec      = elec_aligned;
    cfg.headshape = headmodel.bnd(1);
    elec_aligned  = ft_electroderealign(cfg);
    save(['output/elec_aligned/' pairs{n_sub} '/' subjects{n_sub} '_elec_aligned.mat'],'elec_aligned');
    %{
    1. Set alpha level to 0.9
    2. Goal is to have as many electrodes as possible on the scalp
    3. Take screenshot with transpose values and save it in the eeg_sources
    folder
    %}
    close all
end

%% == 1) Load data and setup parameters ==================================
% Initilize EEGLAB and get participants and pairs names
fileID = fopen('data/pair_codes.csv');
participant_info = textscan(fileID,'%s %s', 'Delimiter', ','); % code, pair
subjects = participant_info{1};
pairs = participant_info{2};
fclose(fileID);
load('dependencies/64ch_withoutA1A2.mat');
ch_names = squeeze(Mon.electrodename);

% Get the duration and time csv files
fileID = fopen('dependencies/play_start_duration.csv');
start_duration = textscan(fileID,'%s %f %s %f %f', 'Delimiter', ',');
fclose(fileID);

% Prepare names of stimuli
conditions = {'hmel', 'hval', 'pclo', 'pcan'};

for n_sub=5:6 %length(subjects)
    mega_eeg = {};
    mega_bad_channels = {};
    mega_conditions = {};
    n_mega = 1;
    for n_condition=1:length(conditions)
        for n_trial=1:5
            % Get file names
            filepath = ['data/raw/' pairs{n_sub} '/eeg/'];
            output = 'output/eeg_eeglab_fieldtrip/';
            names = dir([filepath 'hM_' subjects{n_sub} '_' conditions{n_condition} '*.hdf5']);
            names = {names.name};
            
            for n_filename = 1:length(names)
                if names{n_filename}(16) == num2str(n_trial)
                    filename = names{n_filename};
                    break
                end
            end
            
            chan_locs = [filepath 'hM_EEGDIGI_' subjects{n_sub} '.sfp'];
            
            % Load file
            EEG = pop_loadhdf5('filename', [filepath filename], 'rejectchans', [], 'ref_ch', []);
            for x = 1:(length(EEG.event)-1)
                EEG = pop_editeventvals(EEG,'delete',1);
            end
            
            %% == 2) High pass filter @ 0.5 Hz ==================================
            % default filtering strategy using Hamming window
            EEG = pop_eegfiltnew(EEG, 0.5, []);
            
            %% == 3) Load channel location file ==============================
            EEG = pop_chanedit(EEG, 'load',{chan_locs 'filetype' 'sfp'}, ...
                'changefield',{4 'datachan' 0});
            EEG = pop_select( EEG,'nochannel',{'A-63' 'A-64'});
            chan_labels = EEG.chanlocs;
            
            %% == 4) Trim data from third event to end of data =============
            % Search for the value in the play_start_duration file
            for n_duration = 1:length(start_duration{1})
                if strcmp(start_duration{1}{n_duration}, pairs{n_sub}) & ...
                        strcmp(start_duration{3}{n_duration}, conditions{n_condition}) & ...
                        start_duration{2}(n_duration) == n_trial
                    start = start_duration{4}(n_duration);
                    duration = start_duration{5}(n_duration);
                    EEG = pop_rmdat( EEG, {'Trigger 1'}, ...
                        [(start-3) (start + duration + 3)] ,0); % add 3s just in case
                    break
                end
            end
            
            %% == 5) Common Average References ===========================
            EEG = fullRankAveRef(EEG);
            
            %% == 6) Take out line noise using spectral regression =======
            % I use step size == window size (Makoto suggestion)
            EEG = pop_cleanline(EEG, 'Bandwidth',2,'ChanCompIndices',[1:62], ...
                'ComputeSpectralPower', true, 'LineFrequencies',[60 120], ...
                'NormalizeSpectrum',false,'LineAlpha',0.01, ...
                'PaddingFactor',2,'PlotFigures',false,'ScanForLines',true, ...
                'SignalType','Channels','SmoothingFactor',100,'VerboseOutput',1, ...
                'SlidingWinLength', 4,'SlidingWinStep', 4);
            close all;
            
            %% == 7) Run Artifact Subspace Reconstruction ================
            % Run ASR, use 8 SD because data is very artifactual
            % -> high-pass filtering disabled (we are already doing that)
            % -> noise based rejection disabled
            % -> bad window disabled because we cannot afford to loose data
            EEG = clean_rawdata(EEG, 5, 'off', 0.85, 'off', 8, 'off');
            
            % Get info regarding rejected channels
            bad_channels = cell(size(chan_labels));
            n_channels = size(chan_labels)
            for x = 1:n_channels(2)
                if ~EEG.etc.clean_channel_mask(x)
                    bad_channels{x} = chan_labels(x).labels;
                end
            end
            bad_channels = bad_channels(~cellfun('isempty',bad_channels));
            
            % Write text file with session bad channels
            fid = fopen( 'output/bad_channels.txt', 'a' );
            fprintf(fid, [filename, ',']);
            for x=1:length(bad_channels)
                if x < length(bad_channels)
                    fprintf(fid, '%s,', bad_channels{x});
                else
                    fprintf(fid, '%s\n', bad_channels{x});
                end
            end
            fclose(fid);
            
            %% == 8) Output EEG lab file and save it =====================
            EEGout = pop_saveset(EEG, 'filename',[filename(1:16) '.set'],'filepath', [output pairs{n_sub} '/']);
            mega_bad_channels{n_mega} = bad_channels;
            mega_eeg{n_mega} = EEG;
            mega_condition{n_mega} = filename(4:16);
            n_mega = n_mega + 1;
        end
    end
    %% == 9) Process baseline using the same steps =======================
    names = dir([filepath 'hM_' subjects{n_sub} '_baseline*.hdf5']);
    names = {names.name};
    filename = names{1};
    
    % Load file
    EEG = pop_loadhdf5('filename', [filepath filename], 'rejectchans', [], 'ref_ch', []);
    for x = 1:(length(EEG.event)-1)
        EEG = pop_editeventvals(EEG,'delete',1);
    end
    
    % == High pass filter @ 0.5 Hz ==================================
    % default filtering strategy using Hamming window
    EEG = pop_eegfiltnew(EEG, 0.5, []);
    
    % == Load channel location file ==============================
    EEG=pop_chanedit(EEG, 'load',{chan_locs 'filetype' 'sfp'}, ...
        'changefield',{4 'datachan' 0});
    EEG = pop_select( EEG,'nochannel',{'A-63' 'A-64'});
    chan_labels = EEG.chanlocs;
    
    % == Common Average References ===========================
    EEG = fullRankAveRef(EEG);
    
    % == Take out line noise using spectral regression ===========
    % I use step size == window size (Makoto suggestion)
    EEG = pop_cleanline(EEG, 'Bandwidth',2,'ChanCompIndices',[1:62], ...
        'ComputeSpectralPower', true, 'LineFrequencies',[60 120], ...
        'NormalizeSpectrum',false,'LineAlpha',0.01, ...
        'PaddingFactor',2,'PlotFigures',false,'ScanForLines',true, ...
        'SignalType','Channels','SmoothingFactor',100,'VerboseOutput',1, ...
        'SlidingWinLength', 4,'SlidingWinStep', 4);
    close all;
    
    % == Run Artifact Subspace Reconstruction ========================
    % Run ASR, use 8 SD because data is very artifactual
    % -> high-pass filtering disabled (we are already doing that)
    % -> noise based rejection disabled
    % -> bad window disabled because we cannot afford to loose data
    EEG = clean_rawdata(EEG, 5, 'off', 0.85, 'off', 8, 'off');
    
    % Get info regarding rejected channels
    bad_channels = cell(size(chan_labels));
    n_channels = size(chan_labels)
    for x = 1:n_channels(2)
        if ~EEG.etc.clean_channel_mask(x)
            bad_channels{x} = chan_labels(x).labels;
        end
    end
    bad_channels = bad_channels(~cellfun('isempty',bad_channels));
    
    % Write text file with session bad channels
    fid = fopen( 'output/bad_channels.txt', 'a' );
    fprintf(fid, [filename, ',']);
    for x=1:length(bad_channels)
        if x < length(bad_channels)
            fprintf(fid, '%s,', bad_channels{x});
        else
            fprintf(fid, '%s\n', bad_channels{x});
        end
    end
    fclose(fid);
    
    % == Output EEG lab file and save it =====================
    EEGout = pop_saveset(EEG, 'filename',[filename(1:14) '.set'],'filepath', [output pairs{n_sub} '/']);
    mega_bad_channels{n_mega} = bad_channels;
    mega_eeg{n_mega} = EEG;
    mega_condition{n_mega} = filename(4:14);
    
    
    
    %% == 11) Prepare EEG file ==============================
    mega_eeg = {};
    mega_bad_channels = {};
    mega_conditions = {};
    n_mega = 1;
    for n_condition=1:length(conditions)
        for n_trial=1:5
            % Create the MEGA EEG structure
            filepath = ['output/eeg_eeglab_fieldtrip/' pairs{n_sub} '/'];
            names = dir([filepath 'hM_' subjects{n_sub} '_' conditions{n_condition} '*.set']);
            names = {names.name};
            for n_filename = 1:length(names)
                if names{n_filename}(16) == num2str(n_trial)
                    filename = names{n_filename};
                    break
                end
            end
            mega_eeg{n_mega} = pop_loadset('filename', filename, ...
                'filepath', filepath, 'loadmode', 'all');
            
            % Create the MEGA BAD CHANNELS structure
            load('dependencies/chan_labels', 'chan_labels');
            bad_channels = cell(size(chan_labels));
            n_channels = size(chan_labels);
            for x = 1:n_channels(2)
                if ~mega_eeg{n_mega}.etc.clean_channel_mask(x)
                    bad_channels{x} = chan_labels(x).labels;
                end
            end
            bad_channels = bad_channels(~cellfun('isempty',bad_channels));
            mega_bad_channels{n_mega} = bad_channels;
            
            % Create the MEGA CONDITION structure
            mega_condition{n_mega} = filename(4:16);
            n_mega = n_mega + 1;
        end
    end
    
    % Add baseline
    filepath = ['output/eeg_eeglab_fieldtrip/' pairs{n_sub} '/'];
    names = dir([filepath 'hM_' subjects{n_sub} '_baseline*.set']);
    names = {names.name};
    filename = names{1};
    
    mega_eeg{n_mega} = pop_loadset('filename', filename, 'filepath', ...
        filepath, 'loadmode', 'all');
    
    % Create the MEGA BAD CHANNELS structure
    load('dependencies/chan_labels', 'chan_labels');
    bad_channels = cell(size(chan_labels));
    n_channels = size(chan_labels);
    for x = 1:n_channels(2)
        if ~mega_eeg{n_mega}.etc.clean_channel_mask(x)
            bad_channels{x} = chan_labels(x).labels;
        end
    end
    bad_channels = bad_channels(~cellfun('isempty',bad_channels));
    mega_bad_channels{n_mega} = bad_channels;
    
    % Create the MEGA CONDITION structure
    mega_condition{n_mega} = filename(4:14);
    
    % Do a quick round of interpolation
    for n_trial = 1:length(mega_eeg)
        mega_eeg{n_trial} = pop_interp(mega_eeg{n_trial}, chan_labels, 'spherical');
    end
    
    % Create fieldtrip preprocessing data structure from scratch
    % 1. Flag channels that were rejected in all trials
    bad_electrodes = {};
    bad_electrodes = mega_bad_channels{1};
    for n_elec = 1:(length(mega_bad_channels))
        bad_electrodes = intersect(bad_electrodes, mega_bad_channels{n_elec});
    end
    bad_electrodes = unique(bad_electrodes);
    
    
    % 2. Get smallest value duration from mega_eeg
    min_duration = length(mega_eeg{n_trial}.times);
    for n_trial = 1:(length(mega_eeg))
        if min_duration > length(mega_eeg{n_trial}.times)
            min_duration = length(mega_eeg{n_trial}.times);
            min_duration_arg = n_trial;
        end
    end
    
    %% Create the fieldtrip data structure - shortened version for cov
    data = [];
    
    label = {};
    for n_elec = 1:numel(chan_labels)
        label{n_elec} = chan_labels(n_elec).labels;
    end
    data.label = label';
    
    data.fsample = mega_eeg{1}.srate;
    
    trial = {};
    for n_trial = 1:numel(mega_eeg)
        trial{n_trial} = mega_eeg{n_trial}.data(:, 1:min_duration);
    end
    data.trial =  trial;
    time = {};
    for n_trial = 1:numel(mega_eeg)
        time{n_trial} = mega_eeg{min_duration_arg}.times;
    end
    data.time = time;
    
    % Take out channels that were bad in all trials
    bad_elec_struct = {'all'};
    for n_elec = 1:length(bad_electrodes)
        bad_elec_struct = [bad_elec_struct, ['-' bad_electrodes{n_elec}]];
    end
    cfg = [];
    cfg.channel = bad_elec_struct;
    data = ft_preprocessing(cfg, data);
    
    %% Create the fieldtrip data structure - all the trials
    data_full = [];
    
    label = {};
    for n_elec = 1:numel(chan_labels)
        label{n_elec} = chan_labels(n_elec).labels;
    end
    data_full.label = label'; % cell-array containing strings, Nchan*1
    
    data_full.fsample = mega_eeg{1}.srate; % sampling frequency in Hz, single number
    
    trial = {};
    for n_trial = 1:numel(mega_eeg)
        trial{n_trial} = mega_eeg{n_trial}.data(:, :);
    end
    data_full.trial =  trial;
    time = {};
    for n_trial = 1:numel(mega_eeg)
        time{n_trial} = mega_eeg{n_trial}.times;
    end
    data_full.time = time;
    
    % Take out channels that were bad in all trials
    bad_elec_struct = {'all'};
    for n_elec = 1:length(bad_electrodes)
        bad_elec_struct = [bad_elec_struct, ['-' bad_electrodes{n_elec}]];
    end
    
    cfg = [];
    cfg.channel = bad_elec_struct;
    data_full = ft_preprocessing(cfg, data_full);
    
    elec_to_prune = [bad_elec_struct, {'all', '-LPA', '-NZ', '-RPA', '-Ground', '-A-63', '-A-64'}];
    elec_aligned_pruned = ft_channelselection(elec_to_prune, elec_aligned);
    
    %% == 10) Prepare leadfield model ==============================
    cfg = [];
    cfg.channel = elec_aligned_pruned;
    cfg.headmodel = headmodel;
    cfg.resolution = 1;
    cfg.normalize = 'yes';
    cfg.grid.tight = 'yes';
    cfg.grid.resolution = 1;
    cfg.grid.unit = 'cm';
    cfg.elec = elec_aligned;
    leadfield = ft_prepare_leadfield(cfg);
    
    % Double check the leadfield by plotting it
    figure();
    plot3(leadfield.pos(leadfield.inside,1), leadfield.pos(leadfield.inside,2), ...
        leadfield.pos(leadfield.inside,3), 'k.');
    hold on;
    sens = ft_convert_units(elec_aligned, 'cm');
    ft_plot_sens(sens,'style', 'r*');
    close all
    
    %% == 11) Compute Cortical Patches ==============================
    % Get the AAL atlas and do a corase partition into 19 regions
    atlas_file = 'dependencies/fieldtrip-21092018/template/atlas/aal/ROI_MNI_V4.nii';
    patches = load('dependencies/aal-coarse-19.mat');
    patches = patches.patches;
    
    atlas = ft_read_atlas(atlas_file);
    atlas = ft_convert_units(atlas, leadfield.unit);
    patches_matrix = [];
    patches_matrix.label = {};
    patches_matrix.basis = {};
    patches_matrix.leadfield = {};
    patches_matrix.inside = {};
    patches_matrix.centroid = {};
    
    % Get the basis for each of the patches
    
    for n_patch = 1:length(patches)
        % Select points in anatomical regions that make up the patch
        cfg = [];
        cfg.atlas = atlas;
        cfg.inputcoord = atlas.coordsys;
        cfg.roi = patches(n_patch).patterns;
        mask = ft_volumelookup(cfg, leadfield);
        
        % choose leadfield vertices inside the mask
        inside = leadfield.inside & mask(:);
        
        % Plot this just in case
        figure;
        % plot all inside points
        ft_plot_mesh(leadfield.pos(leadfield.inside,:),'vertexcolor','g');
        hold on;
        % plot patch points
        ft_plot_mesh(leadfield.pos(mask,:));
        close all
        
        % get leadfields in patch
        lf_patch = leadfield.leadfield(inside);
        
        % concatenate intoa single wide matrix
        % n_channels x (3 * points in the patch)
        Hk = [lf_patch{:}];
        
        
        % generate basis for patch
        
        % assume Gamma = I, i.e. white noise
        % otherwise
        %   L = chol(Gamma);
        %   Hk_tilde = L*Hk;
        
        % take SVD of Hk
        % S elements are in decreasing order
        [U,Sk,~] = svd(Hk);
        nsingular = size(Sk,1);
        Uk = [];
        % select the minimum number of singular values
        for j=1:nsingular
            % NOTE if Gamma ~= I
            %   Uk_tilde = L*Uk;
            
            % select j left singular vectors corresponding to largest singular
            % values
            Uk = U(:,1:j);
            
            % compute the representation accuracy
            gammak = trace(Hk'*(Uk*Uk')*Hk)/trace(Hk'*Hk);
            
            % check if we've reached the threshold
            if gammak > 0.85 % this is the eta parameter, or the
                %       representation accuracy, ideally should be close to 1 but it will
                %       also lose its ability to differentiate other patches and resolution
                %       will suffer, see Limpiti2006 for more
                % use the current Uk
                break;
            end
        end
        
        
        % compute centroid
        locs = leadfield.pos(inside,:);
        centroid =  mean(locs,1);
        
        % save all information!
        patches_matrix(n_patch).label = patches(n_patch).name;
        patches_matrix(n_patch).basis = Uk;
        patches_matrix(n_patch).leadfield = Hk;
        patches_matrix(n_patch).inside = inside;
        patches_matrix(n_patch).centroid = centroid;
    end
    
    %% We have a forward model. Now compute the spatial filters
    % Get the covariance matrix
    cov_avg = zeros(length(data.label), length(data.label), length(data.trial));
    
    for n_trial = 1:length(data.trial)
        cov_avg(:, :, n_trial) = data.trial{n_trial} * data.trial{n_trial}';
    end
    
    cov_avg = mean(cov_avg, 3);
    
    % allocate mem
    source = [];
    source.filters = cell(length(patches_matrix), 1);
    source.patch_labels = cell(length(patches_matrix), 1);
    source.patch_centroid = zeros(length(leadfield.inside),3);
    source.inside = false(size(leadfield.inside));
    
    
    % computer filter for each patch
    for n_patch=1:length(patches_matrix)
        fprintf('Computing filter for patch: %s\n',patches_matrix(n_patch).label);
        
        
        % get the patch basis
        Uk = patches_matrix(n_patch).basis;
        
        Yk = Uk'*pinv(cov_avg)*Uk;
        
        
        % Moment orientation is unkown, so we maximize the power
        
        % ordered from smallest to largest
        [V,D] = eig(Yk);
        d = diag(D);
        [~,idx] = sort(d(:),1,'ascend');
        
        % select the eigenvector corresponding to the smallest eigenvalue
        vk = V(:,idx(1));
        
        % compute patch filter weights
        filter = pinv(vk'*Yk*vk)*pinv(cov_avg)*Uk*vk;
        filter = filter';
        
        % save patch filter
        source.filters{n_patch} = filter;
        
        % save patch label for each point
        source.patch_labels{n_patch} = patches_matrix(n_patch).label;
        
        % save centroid
        nverts = sum(patches_matrix(n_patch).inside);
        source.patch_centroid(patches_matrix(n_patch).inside,:) = ...
            repmat(patches_matrix(n_patch).centroid,nverts,1);
        
        
        
        source.inside = leadfield.inside;
        
    end
    
    % Now that you have the patches, we get the data
    % save filters
    leadfield.filter = source.filters;
    leadfield.filter_label = source.patch_labels;
    leadfield.patch_centroid = source.patch_centroid;
    leadfield.inside = source.inside;
    
    % Get source data
    cortical_patches = zeros(length(leadfield.filter), length(data.label));
    for n_patch = 1:length(leadfield.filter)
        cortical_patches(n_patch, :) = leadfield.filter{n_patch};
    end
    
    source_data = [];
    source_data.leadfield = leadfield;
    source_data.trials_labels = mega_condition';
    source_data.bad_electrodes = bad_electrodes;
    source_data.source_label = leadfield.filter_label;
    source_data.sources = cell(length(data.trial), 1);
    source_data.fsample = data_full.fsample;
    source_data.times = data_full.time';
    source_data.chan_label = data_full.label;
    
    for n_trial = 1:length(data_full.trial)
        source_data.sources{n_trial} = cortical_patches * data_full.trial{n_trial};
    end
    
    save(['output/eeg_sources/' pairs{n_sub} '/' subjects{n_sub} '_sources.mat'], 'source_data');
end

%% We first plot the power average on the mri template from each of the cortical patches
for n_patch = 1:length(leadfield.filter)
    filt_seed = leadfield.filter{n_patch};
    
    % create source struct
    source = [];
    source.dim = leadfield.dim;
    source.pos = leadfield.pos;
    source.inside = leadfield.inside;
    source.method = 'average';
    source.avg.pow = zeros(size(leadfield.inside));
    
    for i=1:length(leadfield.inside)
        if leadfield.inside(i)
            source.avg.pow(i) = norm(filt_seed*leadfield.leadfield{i},'fro');
        end
    end
    
    cfgin = [];
    resliced = ft_volumereslice(cfgin, mri);
    cfgin = [];
    cfgin.parameter = 'pow';
    interp = ft_sourceinterpolate(cfgin, source, resliced);
    % source plot
    cfgplot = [];
    cfgplot.method = 'slice';
    cfgplot.funparameter = 'pow';
    cfgplot.maskparameter = cfgplot.funparameter;
    cfgplot.funcolorlim   = [0.4 0.8];
    cfgplot.opacitylim    = [0.4 0.8];
    cfgplot.opacitymap    = 'rampup';
    cfgplot.title = leadfield.filter_label{n_patch};
    ft_sourceplot(cfgplot, interp);
end


%% Patch resolution on anatomical image
% You can also plot the patch resolution on an anatomical image
% equation 14 of Limpiti2006
%function plot_patch_resolution(obj,seed,varargin)
%PLOT_PATCH_RESOLUTION plots patch resolution on anatomical image
%   PLOT_PATCH_RESOLUTION(obj, seed, ...) plots patch
%   resolution on anatomical images. The resolution is computed
%   wrt a seed location and is represented by eq (14) of
%   Limpiti2006
%

% reslice
cfgin = [];
resliced = ft_volumereslice(cfgin, mri);

% create source struct
for n_patch = 1:length(patches_matrix)
    source = [];
    source.dim = leadfield.dim;
    source.pos = leadfield.pos;
    source.inside = leadfield.inside;
    source.method = 'average';
    source.avg.pow = zeros(size(leadfield.inside));
    
    H = patches_matrix(n_patch).leadfield;
    U_seed = patches_matrix(n_patch).basis;
    delta = trace(H'*(U_seed*U_seed')*H)/trace(H'*H);
    source.avg.pow(patches_matrix(n_patch).inside) = delta;
    
    % interpolate
    cfgin = [];
    cfgin.parameter = 'pow';
    interp = ft_sourceinterpolate(cfgin, source, resliced);
    
    % source plot
    cfgplot = [];
    cfgplot.maskparameter = 'mask';
    interp.mask = interp.pow > max(interp.pow(:))*0.5;
    cfgplot.method = 'slice';
    cfgplot.funparameter = 'pow';
    
    ft_sourceplot(cfgplot, interp);
    saveas(gcf,['figs/' patches_matrix(n_patch).label '.png']);
    close all
end