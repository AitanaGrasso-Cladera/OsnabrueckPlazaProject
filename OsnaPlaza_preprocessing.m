%% Automatic cleaning and preprocessing
% This script was developed by Vincent Schmidt.
% Small adjustements, saving the cleaned data and a 90Hz line noise filter
% were added by Debora Nolte.
% Adjusted for Osnabrück Plaza Project by Aitana Grasso Cladera
%% History
% 04.06.25 - AGC: The script now loads the triggers from the external
% hardware system, adds them to the trigger file for fixations/saccades and
% then removes non experimental portions of data (i.e., all datapoints
% around hardware triggers).
%                 Different values for the parameter "BurstCriterion" in
% clean_rawdata were tested, and we decided to maintain 20.
%                 Information was addedd to every EEG file saved in every
% step in order to have a clearer overwiew of the process and files.
% 13.06.25 - AGC: Adjustment of Zapline parameters, line noise was not
% being removed in most of the subjects when using 'line'. Frequency range
% was specified.
%% TO DO
%%
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc.
%% What type of ERP you want to compute?
fixationOnset = 1;
saccadeOnset = 0;
%% Path definition
desktop = 0;
laptop = 1;

if desktop == 1
    dataFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza/Data';
    triggerFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza/triggerFiles';
    hardwareFolder = [];
    addpath('/media/agrassoclade/Aitanas_SSD/Matlab_Toolboxes/eeglab2025.0.0');
    workingFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza';
elseif laptop == 1
    dataFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza/Data';
    triggerFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza/triggerFiles';
    hardwareFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza/hardwareTriggerFile';
    addpath('/Volumes/Aitanas_SSD/Matlab_Toolboxes/eeglab2025.0.0');
    workingFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza';
end
%%
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.'|| ~contains(tmp(pId).name,'xdf')
        continue
    else
        participants(inx).name = tmp(pId).name;
        participants(inx).folder = tmp(pId).folder;
        participants(inx).date = tmp(pId).date;
        inx = inx+1;
    end
end

clear tmp
%%
eeglab
for sub = 1:8%length(participants)
    participantFolder = [dataFolder,filesep,'preproc_',participants(sub).name(1:end-4)];

    if ~exist(participantFolder,'dir')
        mkdir(participantFolder)
    end
    %% load the .xdf file
    excludeMrkrStrms = {'RESTApi-Events - PRO_100'}; % All that is not eeg
    EEG = eeg_load_xdf([dataFolder,filesep,participants(sub).name],'exclude_markerstreams',excludeMrkrStrms);
    %% Electrode renaming and rejection of empty channels 
    % Adjusting the channel locations for the specific setup we used for the
    % recording:
    EEG = pop_chanedit(EEG, 'lookup',[workingFolder,filesep,'Smartfones_EEGlab_polar.elp']);
    % Cleaning channels that do not contain EEG data by default
    alldel = {'AccX','AccY','AccZ','GyroX','GyroY','GyroZ','QuatW','QuatX','QuatY','QuatZ'};
    EEG = pop_select(EEG, 'nochannel', alldel);
 
    % Save the data: always add an 'a' behind the number of automated
    EEG.setname = sprintf('0a_rawChanNames_%s',participants(sub).name(1:end-4));
    EEG.etc.step = 'Step 1: Electrode renaming and rejection of empty channels';
    EEG = eeg_checkset(EEG);
    EEG = pop_editset(EEG, 'setname', sprintf('0a_rawChanNames_%s',participants(sub).name(1:end-4))); 
    EEG = pop_saveset(EEG, 'filename',sprintf('0a_rawChanNames_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    %% Import the trigger file, arrange, and remove non empirical pieces of
    % data
    EEG = pop_loadset(sprintf('0a_rawChanNames_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    hardwareTrigger = readtable([hardwareFolder,filesep,'hardwareTriggers_',participants(sub).name(1:end-4),'.csv'], 'TextType', 'string');
    % Arrange the trigger file
    if fixationOnset
        data = readtable([triggerFolder,filesep,'rawTriggerFileFixations_',participants(sub).name(1:end-4),'.csv'], 'TextType', 'string'); % Assuming the file has headers
        % Remove if there is any nan event
        nonValidTrial = find(isnan(data.givenDuration));
        if ~isempty(nonValidTrial)
            data(nonValidTrial,:) = [];
        end
        % Incorporate hardware triggers into the fixation trigger file
        newTabCol = nan(height(hardwareTrigger), 1);
        hardwareTrigger.("givenDuration") = newTabCol;
        hardwareTrigger.("fixationStart") = newTabCol;
        data = [data;hardwareTrigger];
        data = sortrows(data,"latency");
        % Save the clean trigger file
        writetable(data,[triggerFolder,filesep,'clean_',participants(sub).name(1:end-4),'_fixationTriggers.csv']);
        % Load the trigger file
        EEG = pop_importevent(EEG,'event',fullfile(triggerFolder,['clean_',participants(sub).name(1:end-4),'_fixationTriggers.csv']),'skipline', 1,'append','no');

        for i = 1:length(EEG.event)
            EEG.event(i).latency = EEG.event(i).var1;
            EEG.event(i).type = EEG.event(i).var2;
            EEG.event(i).givenDurations = EEG.event(i).var3;
            EEG.event(i).fixationStart = EEG.event(i).var4;
        end

        EEG.event = rmfield(EEG.event, {'var1', 'var2','var3','var4'});

    else
        data = readtable([triggerFolder,filesep,'rawTriggerFileSaccades_',participants(sub).name(1:end-4),'.csv'], 'TextType', 'string'); % Assuming the file has headers
        % Remove if there is any nan event
        nonValidTrial = find(isnan(data.givenDuration));
        if ~isempty(nonValidTrial)
            data(nonValidTrial,:) = [];
        end
        % Incorporate hardware triggers into the fixation trigger file
        newTabCol = nan(height(hardwareTrigger), 1);
        hardwareTrigger.("givenDuration") = newTabCol;
        hardwareTrigger.("saccadeStart") = newTabCol;
        data = [data;hardwareTrigger];
        data = sortrows(data,"latency");

        data = data(:, {'latency', 'type','saccadeStart','givenDuration'});

        % Save the clean trigger file
        writetable(data,[triggerFolder,filesep,'clean_',participants(sub).name(1:end-4),'_saccadeTriggers.csv']);
        % Load the trigger file
        EEG = pop_importevent(EEG, ...
            'event', fullfile(triggerFolder, ['clean_', participants(sub).name(1:end-4), '_saccadeTriggers.csv']), ...
            'skipline', 1, ...
            'append', 'no');
        for i = 1:length(EEG.event)
            EEG.event(i).latency = EEG.event(i).var1;
            EEG.event(i).type = EEG.event(i).var2;
            EEG.event(i).saccadeStart = EEG.event(i).var3;
            EEG.event(i).givenDurations = EEG.event(i).var4;
        end

        EEG.event = rmfield(EEG.event, {'var1', 'var2','var3','var4'});
    end

    % Remove non experimental parts of data
    eventMark = {EEG.event.type}';
    buffer = (10*EEG.srate)+1;
    % For A
    eTriggerA = max(find(strcmp(eventMark, 'A')));

    EEG = pop_select(EEG,'rmpoint',EEG.times(1):EEG.event(eTriggerA).latency+buffer);
    EEG = eeg_checkset(EEG);

    % For B and C
    eventMark = {EEG.event.type}';
    sTriggerB = min(find(strcmp(eventMark, 'B')));
    eTriggerB = max(find(strcmp(eventMark, 'B')));

    EEG = pop_select(EEG,'rmpoint',EEG.event(sTriggerB).latency-buffer:EEG.event(eTriggerB).latency+buffer);
    EEG = eeg_checkset(EEG);

    eventMark = {EEG.event.type}';
    sTriggerC = min(find(strcmp(eventMark, 'C')));
    eTriggerC = max(find(strcmp(eventMark, 'C')));
    EEG = pop_select(EEG,'rmpoint',EEG.event(sTriggerC).latency-buffer:EEG.event(eTriggerC).latency+buffer);
    EEG = eeg_checkset(EEG);

    % For D
    eventMark = {EEG.event.type}';
    sTriggerD = min(find(strcmp(eventMark, 'D')));
    EEG = pop_select(EEG,'rmpoint',EEG.event(sTriggerD).latency-buffer:EEG.times(end));
    EEG = eeg_checkset(EEG);
    %% Filtering
    % Filter the data parameters adapted from Czeszumski, 2023 (Hyperscanning
    % Maastricht)
    low_pass = 128;
    high_pass = 0.5;
    EEG = pop_eegfiltnew(EEG, high_pass, []); % 0.5 is the lower edge
    EEG = pop_eegfiltnew(EEG, [], low_pass); % 128 is the upper edge

    EEG1 = EEG;

    % Resample (as recommended by Zapline Plus)
    EEG = pop_resample(EEG,500);

    % Remove line noise with zapline
    zaplineConfig=[];
    zaplineConfig.noisefreqs=49.97:.01:50.03; %Alternative: 'line', this was not working for us
    EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG, zaplineConfig); EEG.etc.zapline

    % Save the data: always add ane 'a' behind the number of automated
    EEG = eeg_checkset(EEG);
    if fixationOnset
        EEG.setname = sprintf('1a_triggersFilteredFixation_%s',participants(sub).name(1:end-4));
        EEG.etc.step = 'Step 2: Events added and filtering.';
        EEG = pop_editset(EEG, 'setname', sprintf('1a_triggersFilteredFixation_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('1a_triggersFilteredFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        EEG.setname = sprintf('1a_triggersFilteredSaccade_%s',participants(sub).name(1:end-4));
        EEG.etc.step = 'Step 2: Events added and filtering.';
        EEG = pop_editset(EEG, 'setname', sprintf('1a_triggersFilteredSaccade_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('1a_triggersFilteredSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
    %% Channel removal, data cleaning
    if fixationOnset
        EEG = pop_loadset(sprintf('1a_triggersFilteredFixation_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    else
        EEG = pop_loadset(sprintf('1a_triggersFilteredSaccade_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    end
    full_chanlocs = EEG.chanlocs; % used for data cleaning and interpolation
    EEG.urchanlocs = full_chanlocs;

    % Resample
    EEG = pop_resample(EEG,500);
    % Compute average to mastoids
    EEG = pop_reref( EEG, [4, 8]); % [4, 8] -> L4, R4

    % Clean data using the clean_rawdata plugin
    [EEG,HP,BUR] = clean_artifacts(EEG,'BurstCriterion',20,'BurstRejection','off'); % BurstCriterion at 20 is what we originally used

    % Recompute average reference
    EEG = pop_reref( EEG,[],'interpchan',[]); 
    BUR = pop_reref( BUR,[],'interpchan',[]); 

    Zr=find(EEG.etc.clean_sample_mask == 0); % find all rejected elements
    if ~isempty(Zr)
        starts = Zr(1);
        ends = [];
        for z = 2:length(Zr)
            if Zr(z-1) + 1 ~= Zr(z)
                starts = [starts, Zr(z)];
                ends = [ends, Zr(z-1)];
            end
        end
        ends = [ends, Zr(z)];
        tmprej = [starts;ends]'; % save the noisy segments (beginning & end)
        % save the removed intervals
        %save(fullfile(participantFolder,sprintf('removed_intervals_%s.mat',participants(sub).name(1:end-4))),'tmprej');
    end
    % save removed channels
    removed_channels = ~ismember({full_chanlocs.labels},{EEG.chanlocs.labels});
    EEG.removed_channels = {full_chanlocs(removed_channels).labels};
    % save the removed channels
    save(fullfile(participantFolder,sprintf('removed_channels_%s.mat',participants(sub).name(1:end-4))),'removed_channels');

    % Save the data: always add an 'a' behind the number of automated
    EEG = eeg_checkset(EEG);
    if fixationOnset
        EEG.setname = sprintf('2a_cleanDataChannelsFixation_%s',participants(sub).name(1:end-4));
        EEG.etc.step = 'Step 3: Channel removal, data cleaning.';
        EEG = pop_editset(EEG, 'setname', sprintf('2a_cleanDataChannelsFixation_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('2a_cleanDataChannelsFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        EEG.setname = sprintf('2a_cleanDataChannelsSaccade_%s',participants(sub).name(1:end-4));
        EEG.etc.step = 'Step 3: Channel removal, data cleaning.';
        EEG = pop_editset(EEG, 'setname', sprintf('2a_cleanDataChannelsSaccade_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('2a_cleanDataChannelsSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
    % save the data without the time interval rejection (for Unfold)
    BUR = eeg_checkset(BUR);
    if fixationOnset
        BUR.setname = sprintf('2a_cleanDataChannels_noRejectionFixation_%s',participants(sub).name(1:end-4));
        BUR.etc.step = 'Step 3: Channel removal, data cleaning.';
        BUR = pop_editset(BUR, 'setname', sprintf('2a_cleanDataChannels_noRejectionFixation_%s',participants(sub).name(1:end-4)));
        BUR = pop_saveset(BUR, 'filename',sprintf('2a_cleanDataChannels_noRejectionFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        BUR.setname = sprintf('2a_cleanDataChannels_noRejectionSaccade_%s',participants(sub).name(1:end-4));
        BUR.etc.step = 'Step 3: Channel removal, data cleaning.';
        BUR = pop_editset(BUR, 'setname', sprintf('2a_cleanDataChannels_noRejectionSaccade_%s',participants(sub).name(1:end-4)));
        BUR = pop_saveset(BUR, 'filename',sprintf('2a_cleanDataChannels_noRejectionSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
    
    %% Interpolating missing channels
    if fixationOnset
        EEG = pop_loadset(sprintf('2a_cleanDataChannelsFixation_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
        % Get all channels that need to be interpolated
        EEGchan = pop_loadset(sprintf('1a_triggersFilteredFixation_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    else
        EEG = pop_loadset(sprintf('2a_cleanDataChannelsSaccade_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
        % Get all channels that need to be interpolated
        EEGchan = pop_loadset(sprintf('1a_triggersFilteredSaccade_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    end
    fullChanlocs = EEGchan.chanlocs; % Used for data cleaning and interpolation
    clear EEGchan
    EEG = pop_interp(EEG,fullChanlocs,'spherical');
    EEG = eeg_checkset(EEG);
    % check if duplicate channel label
    if isfield(EEG.chanlocs, 'labels')
        if length({ EEG.chanlocs.labels}) > length(unique({EEG.chanlocs.labels}))
            disp('Warning: some channels have the same label');
        end
    end
    % Save the results
    EEG = eeg_checkset(EEG);
    if fixationOnset
        EEG.setname = sprintf('3a_interpolationFixation_%s',participants(sub).name(1:end-4));
        EEG.etc.step = 'Step 4: Channel interpolation.';
        EEG = pop_saveset(EEG, 'filename',sprintf('3a_interpolationFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        EEG.setname = sprintf('3a_interpolationSaccade_%s',participants(sub).name(1:end-4));
        EEG.etc.step = 'Step 4: Channel interpolation.';
        EEG = pop_saveset(EEG, 'filename',sprintf('3a_interpolationSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
end