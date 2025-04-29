%% Automatic cleaning and preprocessing
% This script was developed by Vincent Schmidt.
% Small adjustements, saving the cleaned data and a 90Hz line noise filter
% were added by Debora Nolte.
% Adjusted for Osnabrück Plaza Project by Aitana Grasso Cladera
%%
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc.
%% What type of ERP you want to compute?
fixationOnset = 0;
saccadeOnset = 1;
%% Path definition
desktop = 0;
laptop = 1;

if desktop == 1
    dataFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza/Data';
    triggerFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza/triggerFiles';
    addpath('/media/agrassoclade/Aitanas_SSD/Matlab_Toolboxes/eeglab2025.0.0');
    workingFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza';
elseif laptop == 1
    dataFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza/Data';
    triggerFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza/triggerFiles';
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
    EEG = eeg_checkset(EEG);
    EEG = pop_editset(EEG, 'setname', sprintf('0a_rawChanNames_%s',participants(sub).name(1:end-4))); 
    EEG = pop_saveset(EEG, 'filename',sprintf('0a_rawChanNames_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    %% Import the trigger file, arrange, filtering 
    EEG = pop_loadset(sprintf('0a_rawChanNames_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    % Arrange the trigger file
    if fixationOnset
        data = readtable([triggerFolder,filesep,'rawTriggerFileFixations_',participants(sub).name(1:end-4),'.csv'], 'TextType', 'string'); % Assuming the file has headers
        % Remove if there is any nan event
        nonValidTrial = find(isnan(data.givenDuration));
        if ~isempty(nonValidTrial)
            data(nonValidTrial,:) = [];
        end
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
        % Save the clean trigger file
        writetable(data,[triggerFolder,filesep,'clean_',participants(sub).name(1:end-4),'_saccadeTriggers.csv']);
        % Load the trigger file
        EEG = pop_importevent(EEG,'event',fullfile(triggerFolder,['clean_',participants(sub).name(1:end-4),'_saccadeTriggers.csv']),'skipline', 1,'append','no');
        for i = 1:length(EEG.event)
            EEG.event(i).latency = EEG.event(i).var2;
            EEG.event(i).type = EEG.event(i).var3;
            EEG.event(i).saccadeStart = EEG.event(i).var1;
            EEG.event(i).givanDurations = EEG.event(i).var4;
        end

        EEG.event = rmfield(EEG.event, {'var1', 'var2','var3','var4'});
    end
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
    zaplineConfig.noisefreqs='line';%49.97:.01:50.03; %Alternative: 'line'
    EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG, zaplineConfig); EEG.etc.zapline

    % Save the data: always add an 'a' behind the number of automated
    EEG = eeg_checkset(EEG);
    if fixationOnset
        EEG = pop_editset(EEG, 'setname', sprintf('1a_triggersFilteredFixation_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('1a_triggersFilteredFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
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
    [EEG,HP,BUR] = clean_artifacts(EEG,'BurstCriterion',20,'BurstRejection','off');

    % [EEG,HP,BUR] = clean_artifacts(EEG,'FlatlineCriterion','off','ChannelCriterion','off',...
    %             'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
    %             'WindowCriterion',0.25,'BurstRejection','off','Distance','Euclidian',...
    %             'WindowCriterionTolerances',[-Inf 7]);
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
        save(fullfile(participantFolder,sprintf('removed_intervals_%s.mat',participants(sub).name(1:end-4))),'tmprej');
    end
    % save removed channels
    removed_channels = ~ismember({full_chanlocs.labels},{EEG.chanlocs.labels});
    EEG.removed_channels = {full_chanlocs(removed_channels).labels};
    % save the removed channels
    save(fullfile(participantFolder,sprintf('removed_channels_%s.mat',participants(sub).name(1:end-4))),'removed_channels');

    % Save the data: always add an 'a' behind the number of automated
    EEG = eeg_checkset(EEG);
    if fixationOnset
        EEG = pop_editset(EEG, 'setname', sprintf('2a_cleanDataChannelsFixation_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('2a_cleanDataChannelsFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        EEG = pop_editset(EEG, 'setname', sprintf('2a_cleanDataChannelsSaccade_%s',participants(sub).name(1:end-4)));
        EEG = pop_saveset(EEG, 'filename',sprintf('2a_cleanDataChannelsSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
    % save the data without the time interval rejection (for Unfold)
    BUR = eeg_checkset(BUR);
    if fixationOnset
        BUR = pop_editset(BUR, 'setname', sprintf('2a_cleanDataChannels_noRejectionFixation_%s',participants(sub).name(1:end-4)));
        BUR = pop_saveset(BUR, 'filename',sprintf('2a_cleanDataChannels_noRejectionFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        BUR = pop_editset(BUR, 'setname', sprintf('2a_cleanDataChannels_noRejectionSaccade_%s',participants(sub).name(1:end-4)));
        BUR = pop_saveset(BUR, 'filename',sprintf('2a_cleanDataChannels_noRejectionSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
    %% ICA - WE ARE SKIPPING THIS BECAUSE IT IS NOT AFFECTING THE DATA - CHECK FOR LATER
    % EEG = pop_loadset(sprintf('2a_cleanDataChannels_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    % 
    % % Create a folder to save the ICA outputs
    % mkdir(fullfile(participantFolder,['amica_',participants(sub).name(1:end-4)]))
    % addpath(fullfile(participantFolder,['amica_',participants(sub).name(1:end-4)]))
    % permission_cleanup(participantFolder);
    % outDir = fullfile(participantFolder,['amica_',participants(sub).name(1:end-4)]);
    % cd(outDir)
    % 
    % 
    % % highpass-filter the data at 2 Hz to not include slow drifts in the ICA
    % eeg_tmp = pop_eegfiltnew(EEG, 2, []);   
    % dataRank = rank(double(eeg_tmp.data'));
    % 
    %  runamica15(eeg_tmp.data, 'num_chans', eeg_tmp.nbchan,'outdir', outDir,... 
    %      'numprocs', 1,'numprocs',1,'max_threads',8,'pcakeep',dataRank,'num_models',1);
    % mod = loadmodout15(ICAfolder);
    % Apply ICA weights to data
    % EEGOUT.icasphere = mod.S;
    % EEGOUT.icaweights = mod.W;
    % EEGOUT.icawinv = [];
    % EEGOUT.icaact = [];
    % EEGOUT.icachansind = [];
    % EEGOUT = eeg_checkset(EEGOUT);
    % Use iclabel to determine which ICs to reject
    % EEGOUT = iclabel(EEGOUT);
    % pop_viewprops(EEG, 0)
    % List components that should be rejected
    % components2remove = [];
    % for component = 1:length(EEGOUT.chanlocs)-1
    %     % Muscle
    %     if EEGOUT.etc.ic_classification.ICLabel.classifications(component,2) > .80
    %         components2remove = [components2remove component];
    %     end
    %     % Eye
    %     if EEGOUT.etc.ic_classification.ICLabel.classifications(component,3) > .9
    %         components2remove = [components2remove component];
    %     end
    %     % Heart
    %     if EEGOUT.etc.ic_classification.ICLabel.classifications(component,4) > .9
    %         components2remove = [components2remove component];
    %     end
    %     % Line noise
    %     if EEGOUT.etc.ic_classification.ICLabel.classifications(component,5) > .9
    %         components2remove = [components2remove component];
    %     end
    %     % Channel noise
    %     if EEGOUT.etc.ic_classification.ICLabel.classifications(component,6) > .9
    %         components2remove = [components2remove component];
    %     end
    % end      
    % % Remove components
    % EEGOUT = pop_subcomp(EEGOUT, components2remove, 0);
    % % Save removed components in struct
    % EEGOUT.removedComponents = components2remove;
    % % Save the removed components
    % save(fullfile(saveFolder,sprintf('removed_components_%s.mat',participants(sub).name(1:end-4))),'components2remove');
    % EEGOUT.etc.Step_4 = 'ICA.';
    % EEGOUT = eeg_checkset(EEGOUT);
    % EEGOUT = pop_saveset(EEGOUT, 'filename',sprintf('4a_ICA_%s',participants(sub).name(1:end-4))),'filepath',fullfile(participantFolder));
    %% Step 5: Interpolating missing channels
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
        EEG = pop_saveset(EEG, 'filename',sprintf('5a_interpolationFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        EEG = pop_saveset(EEG, 'filename',sprintf('5a_interpolationSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
end