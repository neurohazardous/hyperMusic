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

This script intakes EEG data from a gTec system and outputs the
preprocessed, decomposed EEG ready for analysis.

I use EEGLAB and plugins (https://sccn.ucsd.edu/eeglab/index.php) to
do all the preprocessing. I created this pipeline mainly based on
Makoto's preprocessing pipeline
(https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline)
%}

%% == hyperMusic: Preprocessing =========================================
%{
1) Headmodel based on based on ICBM 152 atlas
2) Load file and delete extra events
3) High pass filter 1 Hz (Hamming windowed sinc FIR)
4) Load channel locations, delete extra channels, and assign 10 20 name
5) Take out line noise using spectral regression
5. Run artifact subspace reconstruction
6. Get info regarding interpolated channels
7. Interpolate bad channels using spherical method
8. Perform CAR
9. Trim data around trigger event
10. Save output .set file
11. Write bad channels to txt file
%}

%% == 0) Setup path ======================================================
addpath('dependencies/');
addpath('dependencies/eeglab14_1_2b/');
addpath('dependencies/fieldtrip-21092018/')
ft_defaults()
addpath('dependencies/data-headmodel/mni152')
addpath('dependencies/spm12/')

%% == 8) Headmodel based on based on ICBM 152 atlas ==================
% This bit is the MRI step of the pipeline
% load MRI data for ICMB 152
mri = ft_read_mri('dependencies/MNI152_T1_0.5mm.nii');
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
outputfile = 'output/icbm152_mri.mat';
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
cfg.spmversion='spm12'
cfg.output = {'brain','skull','scalp'};

seg = ft_volumesegment(cfg, mri);

outputfile = 'output/icbm152_seg.mat';
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

outputfile = 'output/icbm152_mesh.mat';
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
mesh = mesh.mesh

cfg = [];
cfg.method = 'dipoli';
cfg.spmversion='spm12';
headmodel = ft_prepare_headmodel(cfg,mesh);
headmodel = ft_convert_units(headmodel,'mm');

outputfile = 'output/icbm152_bem.mat';
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

%% == 1) Load data and setup parameters ==================================
% Initilize EEGLAB and get participants and pairs names
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab
close all
fileID = fopen('data/pair_codes.csv');
participant_info = textscan(fileID,'%s %s', 'Delimiter', ','); % code, pair
subjects = participant_info{1}
pairs = participant_info{2}
fclose(fileID);
load('dependencies/64ch_withoutA1A2.mat');
ch_names = squeeze(Mon.electrodename);

% Get the duration and time csv files
fileID = fopen('dependencies/play_start_duration.csv');
start_duration = textscan(fileID,'%s %f %s %f %f', 'Delimiter', ',');
fclose(fileID);

% Prepare names of stimuli
pieces = {'hmel', 'hval', 'pclo', 'pcan'};

for n_sub=5:5 %length(subjects)
    % Process baseline here!
    for n_piece=1:length(pieces)
        for n_trial=1:5
            % Get file names
            filepath = ['data/raw/' pairs{n_sub} '/eeg/'];
            output = 'output/preprocessed';
            names = dir([filepath 'hM_' subjects{n_sub} '_' pieces{n_piece} '*.hdf5']);
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
            
            %% == 3) High pass filter @ 1Hz ==================================
            % default filtering strategy using Hamming window
            EEG = pop_eegfiltnew(EEG, [], 1, 7920, true, [], 1);
            
            %% == 4) Load channel location file ==============================
            EEG=pop_chanedit(EEG, 'load',{chan_locs 'filetype' 'sfp'}, ...
                'changefield',{4 'datachan' 0});
            EEG = pop_select( EEG,'nochannel',{'A-63' 'A-64'});
            for n = 1:length(ch_names)
                EEG.chanlocs(n).labels = ch_names{n};
            end
            
            % Save original EEG before channel rejection/artifact correction
            originalChanLocs = EEG.chanlocs;
            
            %% == 5) Trim data from third event to end of data =============
            for n_duration = 1:length(start_duration{1})
                if start_duration{1}{n_duration} == pairs{n_sub}
                    if start_duration{3}{n_duration} == pieces{n_piece}
                        if start_duration{2}(n_duration) == n_trial
                            start = start_duration{4}(n_duration);
                            duration = start_duration{5}(n_duration);
                            break
                        end
                    end
                end
                
            end
            
            EEG = pop_rmdat( EEG, {'Trigger 1'}, ...
                [(start-3) (duration + start + 3)] ,0);
            
            %% == 5) Take out line noise using spectral regression ===========
            % I use step size == window size (Makoto suggestion)
            EEG = pop_cleanline(EEG, 'Bandwidth',2,'ChanCompIndices',[1:62], ...
                'ComputeSpectralPower', true, 'LineFrequencies',[60 120], ...
                'NormalizeSpectrum',false,'LineAlpha',0.01, ...
                'PaddingFactor',2,'PlotFigures',false,'ScanForLines',true, ...
                'SignalType','Channels','SmoothingFactor',100,'VerboseOutput',1, ...
                'SlidingWinLength', 4,'SlidingWinStep', 4);
            
            %% == 6) Run Artifact Subspace Reconstruction ========================
            % Run ASR, use 8 SD because data is very artifactual
            % -> high-pass filtering disabled (we are already doing that)
            % -> noise based rejection disabled
            % -> bad window disabled because we cannot afford to loose data
            EEG = clean_rawdata(EEG, 5, [-1], 0.85, -1, 8, -1);
            
            % Get info regarding rejected channels
            bad_channels = cell(1,length(EEG.etc.clean_channel_mask));
            for x = 1:length(EEG.etc.clean_channel_mask)
                if ~EEG.etc.clean_channel_mask(x)
                    bad_channels{x} = ch_names{x};
                end
            end
            bad_channels = bad_channels(~cellfun('isempty',bad_channels));
            
            %% == 7) Common Average References ===========================
            EEG = fullRankAveRef(EEG);
            
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
            EEGout = pop_saveset(EEG, 'filename',[filename(1:16) '.set'],'filepath', [output '/eeg']);
            
        end
    end
    %% == 9) Electrode allignment with template ==============================
    % Properties
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
    ft_plot_sens(elec, 'marker', 's', 'color', 'black');
    
    % create a structure similar to a template set of electrodes
    fid = [];
    fid.elecpos       = [nas; lpa; rpa];       % ctf-coordinates of fiducials
    fid.label         = {'NZ','LPA','RPA'};    % same labels as in elec
    fid.unit          = 'mm';                  % same units as mri
    
    % alignment
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
    
    
    % remove bad electrodes and fiducials TO DO!
    for x = 1:1
        fidch = {'all', ['-' obj.fid_nas],['-' obj.fid_lpa],['-' obj.fid_rpa]};
        channels = ft_channelselection(fidch, sens.label);
    end
    
    %% == 10) Prepare leadfield model ==============================
    cfg = [];
    cfg.headmodel = headmodel;
    cfg.resolution = 1;
    cfg.normalize = 'yes';
    cfg.grid.tight = 'yes';
    cfg.grid.resolution = 1;
    cfg.grid.unit = 'cm';
    cfg.elec = elec_aligned
    leadfield = ft_prepare_leadfield(cfg);
    
    % Double check the leadfield by plotting it
    figure();
    plot3(leadfield.pos(leadfield.inside,1), leadfield.pos(leadfield.inside,2), ...
        leadfield.pos(leadfield.inside,3), 'k.');
    hold on;
    sens = ft_convert_units(elec_aligned, 'cm');
    ft_plot_sens(sens,'style', 'r*');
    
    %% == 11) Prepare EEG file ==============================
    eeg = eeglab2fieldtrip(EEG, 'timelockanalysis', 'none');
    % most of the next code expects the data to be an output of
    % ft_timelockanalysis!!!
    % Plot data
    cfg = [];
    cfg.elec = elec_aligned;
    layout = ft_prepare_layout(cfg, []);
    
    % So, you have two options here -- you'll figure it out once you do the
    % second pass
    cfg = [];
    cfg.layout = layout;
    cfg.showlabels = 'yes';
    ft_multiplotER(cfg, eeg)
    
    % or you can use ft_databrwoser = ft_databrower(cfg, data)
    
    %% == 11) Beamfomer ==============================
    
    % For this step you will need timelocked data!!!
    patchmodel_name = 'aal-coarse-13';
    PathModel = {'aal-coarse-13'};
    cov_avg = 'no';
    compute_lcmv_patch_filters = {'mode', 'single', 'fixedori', 1};
    cfg = []
    cfg.rawtrial = 'yes';
    cfg.method = 'lcmv';
    cfg.grid = leadfield;
    cfg.keepmom = 'yes';
    cfg.lambda = .01;
    cfg.elec = elec_aligned;
    cfg.headmodel = headmodel;
    % You have to feed ft_sourceanalysis to the covariance matrices!
    sourceanalysis = ft_sourceanalysis(cfg, covariances_trials)
    
    % Insert here the normalization using the baseline!
    
    % Plot source power and anatomical image; if you have issues here take a
    % look at here: http://www.fieldtriptoolbox.org/tutorial/beamformer
    % reslice
    cfg = [];
    resliced = ft_volumereslice(cfg, mri);
    % source is supposed to be one of the soruces (I think) you can check
    % fieldtrip later!
    source = rmfield(source,'time');
    % interpolate
    cfgin = [];
    cfgin.parameter = 'pow';
    interp = ft_sourceinterpolate(cfgin, source, resliced);
    
    
    % source plot
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.method = 'slice';
    figure;
    ft_sourceplot(cfg, interp);
    
    % We now remove the outliers from this
    source = sourceanalysis;
    pow = source.avg.pow;
    % create an index
    idx = 1:length(pow);
    temp = [pow(:) idx(:)];
    % sort by source power
    sorted = sortrows(temp,-1);
    sorted(isnan(sorted(:,1)),:) = [];
    
    fprintf('found %d: %f\n', fliplr(sorted(1:n,:))');
    
    % get indices of top n sources
    idx_zero = sorted(1:n,2);
    % zero the top most
    pow(idx_zero) = NaN;
    
    % save data
    source.avg.pow = pow;
    
    
    %% == 11) Compute Cortical Patches ==============================
    % Get the AAL atlas and do a corase partition into 13 regions
    atlas_file = 'dependencies/fieldtrip-21092018/template/atlas/aal/ROI_MNI_V4.nii';
    
    
    patches = [];
    k = 1;
    patches(k).name = 'Prefrontal Left';
    patches(k).patterns = {...
        'Frontal_Sup_L',...
        'Frontal_Sup_Orb_L',...
        'Frontal_Mid_L',...
        'Frontal_Mid_Orb_L',...
        'Frontal_Inf_Oper_L',...
        'Frontal_Inf_Tri_L',...
        'Frontal_Inf_Orb_L',...
        'Frontal_Sup_Medial_L',...
        'Frontal_Med_Orb_L',...
        'Rolandic_Oper_L',...
        'Rectus_L',...
        'Olfactory_L'};
    k = k+1;
    
    patches(k).name = 'Prefrontal Right';
    patches(k).patterns = {...
        'Frontal_Sup_R',...
        'Frontal_Sup_Orb_R',...
        'Frontal_Mid_R',...
        'Frontal_Mid_Orb_R',...
        'Frontal_Inf_Oper_R',...
        'Frontal_Inf_Tri_R',...
        'Frontal_Inf_Orb_R',...
        'Frontal_Sup_Medial_R',...
        'Frontal_Med_Orb_R',...
        'Rolandic_Oper_R',...
        'Rectus_R',...
        'Olfactory_R'};
    k = k+1;
    
    
    patches(k).name = 'Motor Left';
    patches(k).patterns = {...
        'Precentral_L',...
        'Supp_Motor_Area_L'};
    k = k+1;
    
    
    
    patches(k).name = 'Motor Right';
    patches(k).patterns = {...
        'Precentral_R',...
        'Supp_Motor_Area_R'};
    k = k+1;
    
    patches(k).name = 'Basal Ganglia Left';
    patches(k).patterns = {...
        'Pallidum_L',...
        'Caudate_L',...
        'Putamen_L',...
        'Thalamus_L'};
    k = k+1;
    
    patches(k).name = 'Basal Ganglia Right';
    patches(k).patterns = {...
        'Pallidum_R',...
        'Caudate_R',...
        'Putamen_R',...
        'Thalamus_R'};
    k = k+1;
    
    
    
    patches(k).name = 'Insula Left';
    patches(k).patterns = {...
        'Insula_L'};
    k = k + 1;
    
    
    patches(k).name = 'Insula Right';
    patches(k).patterns =  {...
        'Insula_R'};
    k = k+1;
    
    
    patches(k).name = 'Parietal Left';
    patches(k).patterns = {...
        'Parietal_Sup_L',...
        'Parietal_Inf_L',...
        'SupraMarginal_L',...
        'Angular_L',...
        'Precuneus_L',...
        'Postcentral_L',...
        'Paracentral_Lobule_L'};
    k = k+1;
    
    patches(k).name = 'Parietal Right';
    patches(k).patterns = {...
        'Parietal_Sup_R',...
        'Parietal_Inf_R',...
        'SupraMarginal_R',...
        'Angular_R',...
        'Precuneus_R',...
        'Postcentral_R',...
        'Paracentral_Lobule_R'};
    k = k+1;
    
    
    
    patches(k).name = 'Temporal Left';
    patches(k).patterns = {...
        'Temporal_Mid_L',...
        'Temporal_Inf_L',...
        'Temporal_Pole_Sup_L',...
        'Temporal_Pole_Mid_L',...
        'Temporal_Sup_L',...
        'Heschl L'};
    k = k + 1;
    
    
    patches(k).name = 'Temporal Right';
    patches(k).patterns = {...
        'Temporal_Mid_R',...
        'Temporal_Inf_R',...
        'Temporal_Pole_Sup_R',...
        'Temporal_Pole_Mid_R',...
        'Temporal_Sup_R',...
        'Heschl R'};
    k = k+1;
    
    
    patches(k).name = 'Occipital Left';
    patches(k).patterns = {...
        'Occipital_Sup_L',...
        'Occipital_Mid_L',...
        'Occipital_Inf_L',...
        'Cuneus_L',...
        'Fusiform_L',...
        'Lingual_L',...
        'Calcarine_L'};
    k = k+1;
    
    patches(k).name = 'Occipital Right';
    patches(k).patterns = {...
        'Occipital_Sup_R',...
        'Occipital_Mid_R',...
        'Occipital_Inf_R',...
        'Cuneus_R',...
        'Fusiform_R',...
        'Lingual_R',...
        'Calcarine_R'};
    k = k+1;
    
    
    patches(k).name = 'Limbic Left';
    patches(k).patterns = {...
        'Hippocampus_L',...
        'ParaHippocampal_L',...
        'Amygdala_L',...
        'Cingulum_Ant_L',...
        'Cingulum_Mid_L',...
        'Cingulum_Post_L'};
    k = k+1;
    
    
    patches(k).name = 'Limbic Right';
    patches(k).patterns = {...
        'Hippocampus_R',...
        'ParaHippocampal_R',...
        'Amygdala_R',...
        'Cingulum_Ant_R',...
        'Cingulum_Mid_R',...
        'Cingulum_Post_R'};
    k = k+1;
    
    
    patches(k).name = 'Cerebellum Left';
    patches(k).patterns = {...
        'Cerebelum_Crus1_L',...
        'Cerebelum_Crus2_L',...
        'Cerebelum_3_L',...
        'Cerebelum_4 5_L',...
        'Cerebelum_6_L',...
        'Cerebelum_7b_L',...
        'Cerebelum_8_L',...
        'Cerebelum_9_L',...
        'Cerebelum_10_L'};
    k = k+1;
    
    
    patches(k).name = 'Cerebellum Right';
    patches(k).patterns = {...
        'Cerebelum_Crus1_R',...
        'Cerebelum_Crus2_R',...
        'Cerebelum_3_R',...
        'Cerebelum_4 5_R',...
        'Cerebelum_6_R',...
        'Cerebelum_7b_R',...
        'Cerebelum_8_R',...
        'Cerebelum_9_R',...
        'Cerebelum_10_R'};
    k = k+1;
    
    patches(k).name = 'Cerebellum Mid';
    patches(k).patterns = {...
        'Vermis_1_2',...
        'Vermis_3',...
        'Vermis_4_5',...
        'Vermis_6',...
        'Vermis_7',...
        'Vermis_8',...
        'Vermis_9',...
        'Vermis_10'};
    
    
    atlas = ft_read_atlas(atlas_file);
    atlas = ft_convert_units(atlas, leadfield.unit);
    patches_matrix = [];
    patches_matrix.label = {};
    patches_matrix.basis = {};
    patches_matrix.leadfield = {};
    patches_matrix.inside = {};
    patches_matrix.centroid = {};
    
    % Get the basis for each of the patches, starting here, you will need
    % to loop this for each patch!
    
    for n_patch = 1:length(patches)
        % Select points in anatomical regions that make up the patch
        cfg = [];
        cfg.atlas = atlas;
        cfg.inputcoord = atlas.coordsys;
        cfg.roi = patches(n_patch).patterns;
        mask = ft_volumelookup(cfg, leadfield);
        
        % choose leadfield vertices inside the mask
        inside = leadfield.inside & mask(:);
        
        %     % Plot this just in case
        %     figure;
        %     % plot all inside points
        %     ft_plot_mesh(leadfield.pos(leadfield.inside,:),'vertexcolor','g');
        %     hold on;
        %     % plot patch points
        %     ft_plot_mesh(leadfield.pos(mask,:));
        
        
        % get leadfields in patch
        lf_patch = leadfield.leadfield(inside);
        
        % concatenate intoa single wide matrix
        % n_channels x (3 * points in the patch)
        Hk = [lf_patch{:}];
        
        % assume Gamma = I, i.e. white noise
        % otherwise
        %   L = chol(Gamma);
        %   Hk_tilde = L*Hk;
        
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
    
    % We have a forward model. Now compute the spatial filters
    
    % allocate mem
    source = [];
    source.filters = cell(size(leadfield.leadfield));
    source.patch_labels = cell(size(leadfield.leadfield));
    source.inside = false(size(leadfield.inside));
    source.patch_centroid = zeros(length(leadfield.inside),3);
    
    % NOTE when supplying ft_sourceanalysis with filters, i can only specify
    % one per grid point, not per trial, so this function can only operate on a
    % single covariance
    data = eeg;
    
    % We can only compute filters with one cov, so we need to use single
    % trials covariance.
    
    ndims_cov = length(size(data.cov));
    if ndims_cov == 3
        ntrials = size(data.cov,1);
        if ntrials == 1
            data.cov = squeeze(data.cov);
        end
    else
        ntrials = 1;
    end
    
    % computer filter for each patch
    for n_patch=1:length(patches_matrix)
        fprintf('Computing filter for patch: %s\n',patches_matrix(n_patch).label);
        
        if isempty(patches_matrix(n_patch).basis)
            filter = zeros(1,size(data.cov,1));
        else
            % get the patch basis
            Uk = patches_matrix(n_patch).basis;
            
            Yk = Uk'*pinv(data.cov)*Uk;
        end
        
        % Moment orientation is unkown, so we maximize the power
        
        % ordered from smallest to largest
        [V,D] = eig(Yk);
        d = diag(D);
        [~,idx] = sort(d(:),1,'ascend');
        
        % select the eigenvector corresponding to the smallest eigenvalue
        vk = V(:,idx(1));
        
        % compute patch filter weights
        filter = pinv(vk'*Yk*vk)*pinv(data.cov)*Uk*vk;
        filter = filter';
        
        % set patch filter at each point in patch
        [source.filters{patches_matrix(n_patch).inside}] = deal(filter);
        % Yes, this is redundant, but it keeps everything else in Fieldtrip working
        % as normal
        
        % save patch label for each point
        [source.patch_labels{patches_matrix(n_patch).inside}] = deal(patches_matrix(n_patch).label);
        
        % save centroid
        nverts = sum(patches_matrix(n_patch).inside);
        source.patch_centroid(patches_matrix(n_patch).inside,:) = ...
            repmat(patches_matrix(n_patch).centroid,nverts,1);
        
        
        % set empty filters to zero
        
        % check which filters are empty
        grid_empty = cellfun('isempty',source.filters);
        % select those that are inside only
        grid_empty = grid_empty' & leadfield.inside;
        % create a filter with zeros
        filter = zeros(size(filter));
        [source.filters{grid_empty}] = deal(filter);
        
        % check which labels are empty
        grid_empty = cellfun('isempty',source.patch_labels);
        [source.patch_labels{grid_empty}] = deal('');
        
        source.inside = leadfield.inside;
        
    end
    
    % Now that you have the patches, we get the data
    % save filters
    leadfield.filter = source.filters;
    leadfield.filter_label = source.patch_labels;
    leadfield.patch_centroid = source.patch_centroid;
    leadfield.inside = source.inside;
    
    % You can calculate the output of the filters here using
    % ft_sourceanalysis: Besides the source positions, you may also include previously computed
    %   spatial filters and/or leadfields like this
    %     cfg.grid.filter
    %     cfg.grid.leadfield
end

% Some function to perform some plotting so you can check everything is
% good!
% Now we plot the beamformer pattern
% First, choose the seed
filt_seed = leadfield(1).filter

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
source = rmfield(source,'time');
cfgin = [];
cfgin.parameter = 'pow';
interp = ft_sourceinterpolate(cfgin, source, resliced);
% source plot
cfgplot = [];
if ~isempty(p.Results.options)
    % copy options
    cfgplot = copyfields(p.Results.options, cfgplot, fieldnames(p.Results.options));
end

if isfield(cfgplot,'mask') && ~isempty(p.Results.mask)
    warning('overwriting mask field');
end
switch p.Results.mask
    case 'thresh'
        fprintf('creating mask\n');
        cfgplot.maskparameter = 'mask';
        interp.mask = interp.pow > max(interp.pow(:))*p.Results.thresh;
    case 'none'
        % none
end

cfgplot.method = p.Results.method;
cfgplot.funparameter = 'pow';
ft_sourceplot(cfgplot, interp);


% You can also plot the patch resolution on an anatomical image
% equation 14 of Limpiti2006
%function plot_patch_resolution(obj,seed,varargin)
%PLOT_PATCH_RESOLUTION plots patch resolution on anatomical image
%   PLOT_PATCH_RESOLUTION(obj, seed, ...) plots patch
%   resolution on anatomical images. The resolution is computed
%   wrt a seed location and is represented by eq (14) of
%   Limpiti2006
%
%   Inputs
%   ------
%   seed (string or scalar)
%       seed patch, you can specify an patch label or an index
%
%   Parameters
%   ----------
%   method (default = 'slice')
%       plotting method: slice or ortho
%
%   options (struct)
%       options for ft_sourceplot, see ft_sourceplot

p = inputParser;
p.StructExpand = false;
addRequired(p,'seed');
parse(p,seed);

% load patches
patches = ftb.util.loadvar(obj.patches);

% find filter corresponding to the seed
if ischar(p.Results.seed)
    % find a patch corresponding to the label
    match = lumberjack.strfindlisti({patches.name}, p.Results.seed);
    if ~any(match)
        error('could not find %s', p.Results.seed);
    end
    U_seed = patches(match).U;
elseif isscalar(p.Results.seed)
    % use index
    U_seed = patches(p.Results.see).U;
else
    error(['ftb:' mfilename],...
        'unknown seed type');
end

% load leadfield
leadfield = ftb.util.loadvar(obj.lf.leadfield);

% create source struct
source = [];
source.dim = leadfield.dim;
source.pos = leadfield.pos;
source.inside = leadfield.inside;
source.method = 'average';
source.avg.pow = zeros(size(leadfield.inside));

for i=1:length(patches)
    H = patches(i).H;
    delta = trace(H'*(U_seed*U_seed')*H)/trace(H'*H);
    source.avg.pow(patches(i).inside) = delta;
end

% get MRI object
mriObj = obj.get_dep('ftb.MRI');
% load mri data
mri = ftb.util.loadvar(mriObj.mri_mat);

obj.plot_anatomical_deps(mri,source,varargin{:});

%function plot_anatomical_deps(~,mri,source,varargin)
%PLOT_ANATOMICAL_DEPS plots source power on anatomical image
%   PLOT_ANATOMICAL_DEPS(obj, ['method', value, 'options', value]) plots source
%   power on anatomical images. Method can be 'slice' or 'ortho'.
%
%   Input
%   -----
%   mri (struct)
%       mri data, ftb.MRI.mri_mat
%   source (struct)
%       source analysis data, ftb.Beamformer.sourceanalysis
%
%   Parameters
%   ----------
%   method (default = 'slice')
%       plotting method: slice or ortho
%   options (struct)
%       options for ft_sourceplot, see ft_sourceplot
%   mask (default = 'none')
%       mask for functional data, if using this opt
%       thresh - plots values above a threshold
%       none - no mask
%   thresh (default = 0.5)
%       threshold for mask = 'thresh', calculated as a factor
%       of the maximum power

% parse inputs
p = inputParser;
p.StructExpand = false;
addParameter(p,'method','slice',@(x)any(validatestring(x,{'slice','ortho'})));
addParameter(p,'options',[]);
addParameter(p,'mask','none',@(x)any(validatestring(x,{'thresh','none'})));
addParameter(p,'thresh',0.5,@isnumeric);
parse(p,varargin{:});

% reslice
% TODO save instead of redoing
cfgin = [];
resliced = ft_volumereslice(cfgin, mri);

if isfield(source,'time')
    source = rmfield(source,'time');
end


% interpolate
cfgin = [];
cfgin.parameter = 'pow';
interp = ft_sourceinterpolate(cfgin, source, resliced);

% data transformation
plot_log = false;
if plot_log
    interp.pow = db(interp.pow,'power');
end

% source plot
cfgplot = [];
if ~isempty(p.Results.options)
    % copy options
    cfgplot = copyfields(p.Results.options, cfgplot, fieldnames(p.Results.options));
end

if isfield(cfgplot,'mask') && ~isempty(p.Results.mask)
    warning('overwriting mask field');
end
switch p.Results.mask
    case 'thresh'
        fprintf('creating mask\n');
        cfgplot.maskparameter = 'mask';
        interp.mask = interp.pow > max(interp.pow(:))*p.Results.thresh;
    case 'none'
        % none
end

cfgplot.method = p.Results.method;
cfgplot.funparameter = 'pow';
ft_sourceplot(cfgplot, interp);

