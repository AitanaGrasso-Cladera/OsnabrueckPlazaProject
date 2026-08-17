%% Root to mean ratio
%% Paths
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%% For ERP
% Load data from ERPs after averaging sessions from same subject
dataFolder = [projectFolder,filesep,'Data',filesep,'unfold'];
load(fullfile(dataFolder,filesep,'averagedERP.mat'))
% RMS amplitudes (channel 11)
blinkERP = squeeze(mean(finalBlinkOn(:,11,:),1));
saccadeERP = squeeze(mean(finalSaccade(:,11,:),1));
rms_blink = sqrt(mean(blinkERP.^2));
rms_saccade = sqrt(mean(saccadeERP.^2));
% Compute Ratio
rms_ratioERP = rms_blink / rms_saccade;
%% For Time Frequency
% Load data from TFRs after averaging sessions from same subject
load(fullfile(dataFolder,filesep,'averagedTF.mat'))
% RMS amplitudes (channel 11)
blinkERP = squeeze(mean(finalBlinkOn(:,11,:),1));
saccadeERP = squeeze(mean(finalSaccade(:,11,:),1));
rms_blink = sqrt(mean(blinkERP.^2));
rms_saccade = sqrt(mean(saccadeERP.^2));
% Compute Ratio
rms_ratioPower = rms_blink / rms_saccade;

save([dataFolder,filesep,'root2mean_ratios.mat'],'rms_ratioERP','rms_ratioPower')
