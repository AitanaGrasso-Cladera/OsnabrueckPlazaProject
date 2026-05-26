%% Epoch the data according to different events of interest
% This script was developed by Aitana Grasso Cladera
%% History
% 12.02.26 - AGC = Script was modified, rejection was removed and done in a
% separated script.
%%
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc.
%% Path definition
projectFolder = [];
dataFolder = [projectFolder,filesep,'Data'];
addpath('[]/eeglab2026.0.0');
workingFolder = [];
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
close
%%
for sub = 1:length(participants)

    participantFolder = [dataFolder,filesep,'preproc_',participants(sub).name(1:end-4)];

    if ~exist(participantFolder,'dir')
        mkdir(participantFolder)
    end

    EEG = pop_loadset(sprintf('3a_interpolationALL_%s.set',participants(sub).name(1:end-4)),fullfile(participantFolder));

    % Epoch the data
    epochLength = [-0.5 0.7];

    epochFixation = pop_epoch(EEG,{'fixationOnset'},epochLength,'newname','epoched','epochinfo','yes');
    epochFixation.setname = sprintf('4a_epochedFixation_%s',participants(sub).name(1:end-4));
    epochFixation.etc.step = 'Step 5: Epoched to fixation onset.';
    epochFixation = pop_saveset(epochFixation, 'filename',sprintf('4a_epochedFixation_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));

    epochSaccade = pop_epoch(EEG,{'saccadeOnset'},epochLength,'newname','epoched','epochinfo','yes');
    epochSaccade.setname = sprintf('4a_epochedSaccade_%s',participants(sub).name(1:end-4));
    epochSaccade.etc.step = 'Step 5: Epoched to saccade onset.';
    epochSaccade = pop_saveset(epochSaccade, 'filename',sprintf('4a_epochedSaccade_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));

    epochBlinkOnset = pop_epoch(EEG,{'blinkOnset'},epochLength,'newname','epoched','epochinfo','yes');
    epochBlinkOnset.setname = sprintf('4a_epochedBlink_%s',participants(sub).name(1:end-4));
    epochBlinkOnset.etc.step = 'Step 5: Epoched to blink onset.';
    epochBlinkOnset = pop_saveset(epochBlinkOnset, 'filename',sprintf('4a_epochedBlinkOnset_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));

    epochBlinkOffset = pop_epoch(EEG,{'blinkOffset'},epochLength,'newname','epoched','epochinfo','yes');
    epochBlinkOffset.setname = sprintf('4a_epochedBlink_%s',participants(sub).name(1:end-4));
    epochBlinkOffset.etc.step = 'Step 5: Epoched to blink onset.';
    epochBlinkOffset = pop_saveset(epochBlinkOffset, 'filename',sprintf('4a_epochedBlinkOffset_%s',participants(sub).name(1:end-4)),'filepath',fullfile(participantFolder));
end
%%
rmpath('[]/eeglab2026.0.0');
