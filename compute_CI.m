function [confInterval,ERP] = compute_CI(Data,alpha,channel)
%% Function to compute Confidence Intervals over participants
% Usage: confInterval = compute_CI(Data,0.05,43,1);
%% Input
% Data: 3D matrix with ERP per participant. The dimensions should be:
% participants x channels x timepoints
% Alpha: chosen significance value. Default is 0.05
% Channel: target channels to compute the ERP and its respective CI.
% Default is all channels
%% Output
% confInterval = 2D (or 3D if you have more than one channel) matrix with 
% lower (1,:) and upper (2,:) confidence interval for every time point.
% Electrodes x Times x CI
% ERP = 2D matrix with channels x time
%% Created by Aitana Grasso-Cladera, MSc
% October 2024
% aitana.grasso.cladera@uni-osnarueck.de
%% Checking
if ~exist('Data','var')
    error('You need to input the data.')
elseif ndims(Data) < 3
    error('Your matrix does not match the requirement of this function')
end

if ~exist('alpha','var')
    alpha = 0.05;
end

if ~exist('channel', 'var') || isempty(channel)
    fprintf('Computing CI over all electrodes\n');
    channel = 1:size(Data, 2); % Default to all channels
end
%% Get parameters
nParticipants = size(Data,1);
nChannels = size(Data,2);
nTimePoints = size(Data,3);
%% Compute CI
tCrit = tinv(1-alpha/2,nParticipants-1); % Critical value from T distribution
if nChannels == 1
    confInterval = nan(2,nTimePoints);
    % Calculate mean ERP amplitude for each timepoint accross participants
    meanERP = mean(squeeze(Data(:,:)),1);
    % Calculate the SD for each timepoing accros participants
    STD = std(squeeze(Data(:,:)),0,1);
    % Calculate the Standard Error of the mean for each timepoint
    SE = STD/sqrt(nParticipants);
    % Calculate the 95% confidence interval for each timepoint
    confInterval(1,:) = meanERP - tCrit * SE;
    confInterval(2,:) = meanERP + tCrit * SE;
elseif nChannels > 1
    confInterval = nan(nChannels,nTimePoints,2);
    ERP = nan(nChannels,nTimePoints);
    for ch = 1:nChannels
        % Calculate mean ERP amplitude for each timepoint accross participants
        meanERP = mean(squeeze(Data(:,ch,:)),1);
        % Calculate the SD for each timepoing accros participants
        STD = std(squeeze(Data(:,ch,:)),0,1);
        % Calculate the Standard Error of the mean for each timepoint
        SE = STD/sqrt(nParticipants);
        % Calculate the 95% confidence interval for each timepoint
        confInterval(ch,:,1) = meanERP - tCrit * SE;
        confInterval(ch,:,2) = meanERP + tCrit * SE;

        % Save the ERP
        ERP(ch,:) = meanERP;
        clear meanERP STD SE
    end
end
end