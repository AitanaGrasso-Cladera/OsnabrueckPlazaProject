%% Epoch the data and remove trials
% This script was developed by Aitana Grasso Cladera
%% History
% 04.06.25 - AGC = Following Ladouce et al., 2023, we removed aberrant
% activity (threshold of 5 std around the median voltage activity recorded
% across all epochs) were discarded.
% 26.06.25 - AGC = Removing the same trials for saccade onset than for
% fixation onset was tested. The performance was not good, since the data
% did not improved. Hence, we run the same calculation but for saccade
% onset ERPs, which generates different trials to be removed.
%% To do
% 04.06.25 - AGC = So far, we are epoching all trials, no difference for 
% type of categorie of the fixation (e.g., face, backgorund).
%%
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc.
%% Path definition
desktop = 0;
laptop = 1;

if desktop == 1
    dataFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza/Data';
    addpath('/media/agrassoclade/Aitanas_SSD/Matlab_Toolboxes/eeglab2025.0.0');
    workingFolder = '/media/agrassoclade/Aitanas_SSD/Osnabruck_Plaza';
elseif laptop == 1
    dataFolder = '/Volumes/Aitanas_SSD/Osnabruck_Plaza/Data';
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
%% Type of ERP
fixation = 0;
saccade = 1;
%%
eeglab
for sub = 1:8%length(participants)
    participantFolder = [dataFolder,filesep,'preproc_',participants(sub).name(1:end-4)];

    if ~exist(participantFolder,'dir')
        mkdir(participantFolder)
    end
    
    % Load EEG data
    if fixation == 1
        EEG = pop_loadset(sprintf('3a_interpolationFixation_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    else
        EEG = pop_loadset(sprintf('3a_interpolationSaccade_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));
    end
    % Epoch the data
    epochLength = [-0.3 0.5];
    epochEEG = pop_epoch(EEG,{},epochLength,'newname','epoched','epochinfo','yes');

    if fixation == 1
        epochEEG.setname = sprintf('4a_epochedFixation_%s',participants(sub).name(1:end-4));
        epochEEG.etc.step = 'Step 5: Epoched to fixation onset.';
        epochEEG = pop_saveset(epochEEG, 'filename',sprintf('4a_epochedFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    else
        epochEEG.setname = sprintf('4a_epochedSaccade_%s',participants(sub).name(1:end-4));
        epochEEG.etc.step = 'Step 5: Epoched to saccade onset.';
        epochEEG = pop_saveset(epochEEG, 'filename',sprintf('4a_epochedSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
    end
    %% Remove bad trials (5 std as in Ladouce et al., 2024)
    % Reshape EEG data to 1D vector (pool all voltages from all epochs)
        all_data = epochEEG.data(:);
        data_median = median(all_data);
        data_std = std(all_data);
        % Compute the threshold
        threshold = 5 * data_std;
        % Loop through epochs to find max deviation from the median
        num_epochs = size(epochEEG.data, 3);
        bad_epochs = false(1, num_epochs);

        for i = 1:num_epochs
            epoch_data = epochEEG.data(:,:,i);
            max_deviation = max(abs(epoch_data(:) - data_median));
            if max_deviation > threshold
                bad_epochs(i) = true;
            end
        end
        % Save bad epochs
        if fixation == 1
            save(fullfile(participantFolder,sprintf('removed_epochsFixation_%s.mat',participants(sub).name(1:end-4))),'bad_epochs');
            % Reject bad epochs
            epochEEG = pop_rejepoch(epochEEG,bad_epochs,0);
            epochEEG.setname = sprintf('5a_epochedFixationRejection_%s',participants(sub).name(1:end-4));
            epochEEG.etc.step = 'Step 6: Rejection of bad epochs.';
            epochEEG = pop_saveset(epochEEG, 'filename',sprintf('5a_epochedFixationRejection_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
        else
            save(fullfile(participantFolder,sprintf('removed_epochsSaccade_%s.mat',participants(sub).name(1:end-4))),'bad_epochs');
            epochEEG = pop_rejepoch(epochEEG,bad_epochs,0);
            epochEEG.setname = sprintf('5a_epochedSaccadeRejection_%s',participants(sub).name(1:end-4));
            epochEEG.etc.step = 'Step 6: Rejection of bad epochs.';
            epochEEG = pop_saveset(epochEEG, 'filename',sprintf('5a_epochedSaccadeRejection_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
        end
end

