%% Automatic data aligning and creation of trigger file
% This script was developed by Aitana Grasso-Cladera and Debora Nolte
%% History
% 20.03.25 - AGC = Automatization of the script to run over multiple
% subjects. At the moment, only creates the trigger file for
% fixation onset ERPs.
%% To do
% Compute saccade onset ERPs
% We are not checking for events that might be dropped
%% Path definition
desktop = 0;
laptop = 1;

if desktop == 1
    EEGdataFolder = [];
    fixationFolder = [];
    ETeventsFolder = [];
    triggerFolder = [];
   % addpath('EEGLAB path');
elseif laptop == 1
    EEGdataFolder = [];
    fixationFolder = [];
    ETeventsFolder = [];
    triggerFolder = [];
    % addpath('EEGLAB path');
end
%% Get information about participant's files
tmp = dir(fullfile(EEGdataFolder));
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
for sub = 1:length(participants)
    %% Create participant folder
    participantFolder = [EEGdataFolder,filesep,'preproc_',participants(sub).name(1:end-4)];

    if ~exist(participantFolder,'dir')
        mkdir(participantFolder)
    end
    %% Load the data from ET events
    eventsPL = readtable([ETeventsFolder,filesep,'events_',participants(sub).name(end-7:end-4)], 'TextType', 'string');
    % Load the XDF file from the Smarting PRO App
    streams  = load_xdf([EEGdataFolder,filesep,participants(sub).name]);
    % Find and save the RESTApi streams
    for i = 1:length(streams)
        if contains(streams{1,i}.info.name,'RESTApi')
            restAPIdata = streams{1,i};
        end
    end
    %% Remove hardware triggers from the ET event file (if any)
    hardPositions = contains(eventsPL.name,'hardsync');
    eventsPL(hardPositions,:) = [];
    %% Get start of from both timelines (EEG.start message on ET time)
    EEGstartMessage = find(contains(eventsPL.name,'EEG.start'));
    EEGendMessage = find(contains(eventsPL.name,'EEG.end'));
    startET = eventsPL.timestamp_ns_(EEGstartMessage);
    startEEG = eventsPL.name(EEGstartMessage);
    % Since it is a string, get only the numeric part
    numericPart = regexp(startEEG, '\d+(\.\d+)?', 'match');
    startEEG = str2double(numericPart{1});
    %% Get the start message on the RESTApi streams
    EEGstartMessageRESTA = find(contains(restAPIdata.time_series,'EEG.start'));
    EEGendMessageRESTA = find(contains(restAPIdata.time_series,'EEG.end'));
    startRestaPI = restAPIdata.time_stamps(EEGstartMessageRESTA);
    %% Get the ET events between start and end of EEG recording
    ETtimes = eventsPL.timestamp_ns_(EEGstartMessage:EEGendMessage);
    % Set ET timestamps to seconds
    firstETtimes = ETtimes/1e9;
    % Set the first timestamp to zero
    ETtimeInSeconds = firstETtimes - firstETtimes(1);
    %% Compute difference between the number of EEG events on the ET and
    % RESTApi/EEG data
    differenceNumberMessages = length(ETtimeInSeconds) - length(restAPIdata.time_stamps);
    %% Take the EEG data from the first sync event after start event
    EEGtimes = restAPIdata.time_stamps(EEGstartMessageRESTA+1:end);
    % Set EEG time stamps to zero
    EEGtimeInSeconds = EEGtimes - EEGtimes(1);
    %% Compute residual
    residualsTimeline = nan(size(EEGtimeInSeconds));
    for i = 1:length(EEGtimeInSeconds)
        residualsTimeline(i) =  EEGtimeInSeconds(i) - ETtimeInSeconds(i);
    end
    %% Fit a regression to later compute drift
    % Fit a first-degree polynomial (linear regression)
    p = polyfit(1:length(residualsTimeline),residualsTimeline,1);
    % Generate fitted y-values
    fitted_values = polyval(p,1:length(residualsTimeline));
    %% Compute the linear drift of the fitted line
    linearDrift = fitted_values(end)-fitted_values(1);
    %% Create the trigger file (based on fixations)
    % Load fixation on faces/background (Face Detection Alg)
    fixEvents = readtable([fixationFolder,filesep,'fixations_on_face_',participants(sub).name(end-7:end-4)]);
    % Set ET timestamps from fixations to seconds
    fixEventsStart = fixEvents.startTimestamp_ns_/1e9;
    % Set the first timestamp to zero
    fixEventsTime = fixEventsStart - firstETtimes(1);
    %% Remove the fixation information before the EEG started
    times2keep = find(fixEventsTime >= 0);
    fixationTimes = fixEventsTime(times2keep);
    fixEventsCutted = fixEvents(times2keep,:);
    %% Cut the end of the EEG recording
    times2keepEnd = find(fixationTimes <= EEGtimeInSeconds(end));
    fixationTimes = fixationTimes(times2keepEnd);
    fixEventsCutted = fixEvents(times2keepEnd,:);
    fixationTimes = round(fixationTimes,3);
    %% Create an interporalted timeline
    fs = 1000;
    interpTime = linspace(EEGtimeInSeconds(1),EEGtimeInSeconds(end),(EEGtimeInSeconds(end)-EEGtimeInSeconds(1))*fs+1);
    interpTime = round(interpTime,3); % for precision
    %% Find fixations times on the inteprolated time
    actualFixations = nan(size(fixationTimes));
    for i = 1:length(fixationTimes)
        actualFixations(i) = find(interpTime == fixationTimes(i));
    end
    %% Add the drift factor to every sample
    driftValues = linspace(0,linearDrift,length(interpTime));
    timeline = interpTime + driftValues;
    %% Select timepoints that correspond to fixations
    newFixationsTime = timeline(actualFixations);
    % Now we add the first EEG timestamps, to get the data in EEG latencies
    finalFixationsTime = (newFixationsTime*500)+1;
    %% Organize fixation trigger file
    fixEventsCutted.latency = finalFixationsTime';
    for i = 1:length(fixEventsCutted.latency)
        if contains(fixEventsCutted.fixationOnFace(i),'true')
            fixEventsCutted.type(i) = {'face'};
        else
            fixEventsCutted.type(i) = {'background'};
        end
    end
    writetable(fixEventsCutted,[triggerFolder,filesep,'rawTriggerFile_Participant',participants(sub).name(end-7:end-4),'.csv']);

    save([EEGdataFolder,filesep,'preproc_',participants(sub).name(1:end-4),filesep,'alignmentMetadata',participants(sub).name(1:end-4),'.mat'],'differenceNumberMessages','residualsTimeline','linearDrift')
end
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc.