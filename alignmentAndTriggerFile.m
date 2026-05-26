%% Automatic data aligning and creation of trigger file
% This script was developed by Aitana Grasso-Cladera and Debora Nolte
%% History
% 20.03.25 - AGC = Automatization of the script to run over multiple
% subjects. At the moment, only creates the trigger file for
% fixation onset ERPs.
% 02.04.25 - AGC = Automatizated script computes the tnotrigger file for
% saccade onset for multiple subjects. 
% 08.04.25 - AGC - DN = We checked how the durations and latencies were
% assigned, and modified the script to be more accurate.
% 03.06.25 - DN - AGC = We checked saccade vs fixation onset time limits
% for trigger file.
% 30.01.26 - AGC = Added amplitude to the saccade information. Removed Face 
% and Background elements.
% 14.04.26 - AGC = Included blinks as part of the trigger file creation and
% did small changes for the normalization of the timelines.
%% Path definition
projectFolder = [];
EEGdataFolder = [projectFolder,filesep,'Data'];
fixationFolder = [projectFolder,filesep,'fixationFiles'];
saccadeFolder = [projectFolder,filesep,'saccadeFiles'];
blinkFolder = [projectFolder,filesep,'blinkFiles'];
ETeventsFolder = [projectFolder,filesep,'pupilLabEvents'];
triggerFolder = [projectFolder,filesep,'triggerFiles'];

addpath(genpath('[]/eeglab2026.0.0'));
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
close
%%
for sub = 1:length(participants)
    % Create participant folder
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
    hardPositions = contains(eventsPL.name,'A');
    eventsPL(hardPositions,:) = [];
    hardPositions = contains(eventsPL.name,'B');
    eventsPL(hardPositions,:) = [];
    hardPositions = contains(eventsPL.name,'C');
    eventsPL(hardPositions,:) = [];
    hardPositions = contains(eventsPL.name,'D');
    eventsPL(hardPositions,:) = [];
    %% Get the start message on the RESTApi streams
    EEGstartMessageRESTA = find(contains(restAPIdata.time_series,'EEG.start'));
    EEGendMessageRESTA = find(contains(restAPIdata.time_series,'EEG.end'));
    startRestaPI = restAPIdata.time_stamps(EEGstartMessageRESTA);
    endRestaPI = restAPIdata.time_stamps(end);
    %% Get start of from both timelines (EEG.start message on ET time)
    EEGstartMessage = find(contains(eventsPL.name,'EEG.start'));
    % First EEG sync message
    EEGsyncMessage = find(contains(eventsPL.name, 'EEG.sync'), 1);
    EEGendMessage = find(contains(eventsPL.name,'EEG.end'));
    % Determine how to proceed
    if EEGstartMessage < EEGsyncMessage
        EEGstartMessage = EEGstartMessage+1;
        startET = eventsPL.timestamp_ns_(EEGstartMessage);
        startEEG = eventsPL.name(EEGstartMessage);
        % Since it is a string, get only the numeric part
        numericPart = regexp(startEEG, '\d+(\.\d+)?', 'match');
        startEEG = str2double(numericPart{1});
        startMatch = floor(startRestaPI) == floor(startEEG);
        % Get the ET events between start and end of EEG recording
        ETtimes = eventsPL.timestamp_ns_(EEGstartMessage:EEGendMessage);
    else
        EEGstartMessage = EEGstartMessage-1;
        startET = eventsPL.timestamp_ns_(EEGstartMessage);
        startEEG = eventsPL.name(EEGstartMessage);
        % Delete the start event to not have a shift
        eventsPL(EEGstartMessage,:) = [];
        % Since it is a string, get only the numeric part
        numericPart = regexp(startEEG, '\d+(\.\d+)?', 'match');
        startEEG = str2double(numericPart{1});
        startMatch = floor(startRestaPI) == floor(startEEG);
        % Get the ET events between start and end of EEG recording
        ETtimes = eventsPL.timestamp_ns_(EEGstartMessage:EEGendMessage);
    end
    % Set ET timestamps to seconds
    firstETtimes = ETtimes/1e9;
    % Set the first timestamp to zero
    ETtimeInSeconds = firstETtimes - firstETtimes(1);
    %% Take the EEG data from the first sync event after start event
    EEGtimes = restAPIdata.time_stamps(EEGstartMessageRESTA+1:end);
    % Set EEG time stamps to zero
    EEGtimeInSeconds = EEGtimes - EEGtimes(1);
    %% Compute difference between the number of EEG events on the ET and
    % RESTApi/EEG data
    differenceNumberMessages = length(ETtimeInSeconds) - length(restAPIdata.time_stamps);
    %% Compute residual
    residualsTimeline = nan(size(EEGtimeInSeconds));
    for i = 1:length(EEGtimeInSeconds)
        residualsTimeline(i) =  EEGtimeInSeconds(i) - ETtimeInSeconds(i);
    end
    %% Fit a regression to later compute drift
    % Fit a first-degree polynomial (linear regression)
    p = polyfit(1:length(residualsTimeline),residualsTimeline,1);
    intercept = p(2);
    % Generate fitted y-values
    fitted_values = polyval(p,1:length(residualsTimeline));
    %% Compute the linear drift of the fitted line
    linearDrift = fitted_values(end)-fitted_values(1);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Create the trigger file (based on fixations)
    % Load fixation on faces/background (Face Detection Alg)
    fixEvents = readtable([fixationFolder,filesep,'finalFixations',filesep,'fixations_',participants(sub).name(end-7:end-4)]);
   % Cut the fixation data between EEG start and end
    lowerLimitETfix = min(find(fixEvents.startTimestamp_ns_ >= eventsPL.timestamp_ns_(EEGstartMessage)));
    upperLimitETfix = max(find(fixEvents.startTimestamp_ns_ <= eventsPL.timestamp_ns_(EEGendMessage)));
    targetLower = fixEvents.startTimestamp_ns_(lowerLimitETfix);
    targetUpper = fixEvents.startTimestamp_ns_(upperLimitETfix);
    fixEvents = fixEvents(lowerLimitETfix:upperLimitETfix,:);
   % Set ET timestamps from fixations to seconds
    fixEventsStart = fixEvents.startTimestamp_ns_/1e9;
   % Set the first timestamp to zero
    fixEventsTime = fixEventsStart - firstETtimes(1);
   % Round for precision
    fixationTimes = round(fixEventsTime,3);
    %% Create the trigger file (based on saccades)
    % Load the saccade file from Pupil Cloud
    saccEvents = readtable([saccadeFolder,filesep,'newSaccades',filesep,'saccades_',participants(sub).name(end-7:end-4)]);
    % Cut the saccade data between EEG start and end
    lowerLimitETsacc = find(saccEvents.startTimestamp_ns_ >= targetLower, 1, 'first');
    upperLimitETsacc = find(saccEvents.startTimestamp_ns_ <= targetUpper, 1, 'last');
    saccEvents = saccEvents(lowerLimitETsacc:upperLimitETsacc,:);
    % Set ET timestamps from saccades to seconds
    saccEventsStart = saccEvents.startTimestamp_ns_/1e9;
    % Set the first timestamp to zero
    saccEventsTime = saccEventsStart - firstETtimes(1);
    % Round for precision
    saccadeTimes = round(saccEventsTime,3);
    %% Create an interporalted timeline
    fs = 1000;
    interpTime = linspace(EEGtimeInSeconds(1),EEGtimeInSeconds(end),(EEGtimeInSeconds(end)-EEGtimeInSeconds(1))*fs+1);
    interpTime = round(interpTime,3); % for precision
    %% Find fixations times on the inteprolated time
    actualFixations = nan(size(fixationTimes));
    for i = 1:length(fixationTimes)
        if fixationTimes(i) <= interpTime(end)
            actualFixations(i) = find(interpTime == fixationTimes(i));
        else
            actualFixations(i) = NaN;
        end
    end
    actualFixations(isnan(actualFixations)) = [];
    %% Find saccade times on the interpolated time
    actualSaccades = nan(size(saccadeTimes));
    for i = 1:length(saccadeTimes)
        if saccadeTimes(i) <= interpTime(end)
            if isempty(find(interpTime == saccadeTimes(i)))
                actualSaccades(i) = NaN;
            else
                actualSaccades(i) = find(interpTime == saccadeTimes(i));
            end
        else
            actualSaccades(i) = NaN;
        end
    end
    actualSaccades(isnan(actualSaccades)) = [];
    %% Add the drift factor to every sample
    driftValues = intercept + linspace(0,linearDrift,length(interpTime));
    timeline = interpTime + driftValues;
    %% Select timepoints that correspond to fixations
    newFixationsTime = timeline(actualFixations);
    % Now we add the first EEG timestamps, to get the data in EEG latencies
    finalFixationsTime = (newFixationsTime*500)+1;
    %% Select timepoints that correspond to saccades
    newSaccadesTime = timeline(actualSaccades);
    % Now we add the first EEG timestamps, to get the data in EEG latencies
    finalSaccadesTime = (newSaccadesTime*500)+1;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Organize trigger files
    nRows = length(finalFixationsTime);  % number of rows
    % Preallocate each column
    saccEventsCutted = table(nan(nRows,1),nan(nRows,1),strings(nRows,1), nan(nRows,1), ...              % givenDuration
        'VariableNames', {'start','latency','type','givenDuration'});
    fixEventsCutted = table(nan(nRows,1),nan(nRows,1),strings(nRows,1),nan(nRows,1), ...
        'VariableNames', {'start','latency','type','givenDuration'});

    for i = 1:nRows
        fixEventsCutted.start(i) = newFixationsTime(i);
        fixEventsCutted.latency(i) = finalFixationsTime(i);
        fixEventsCutted.type(i) = fixEvents.type(i);
        fixEventsCutted.givenDuration(i) = fixEvents.duration_ms_(i);

        validIndices = find(newSaccadesTime < newFixationsTime(i));

        if ~isempty(validIndices)
            [~, idx] = min(newFixationsTime(i) - newSaccadesTime(validIndices));
            closestIndex = validIndices(idx);

            saccEventsCutted.start(i) = newSaccadesTime(closestIndex);
            saccEventsCutted.latency(i) = finalSaccadesTime(closestIndex);
            saccEventsCutted.type(i) = saccEvents.type(closestIndex);
            saccEventsCutted.givenDuration(i) = saccEvents.duration_ms_(closestIndex);

        else
            saccEventsCutted.start(i) = NaN;
            saccEventsCutted.latency(i) = NaN;
            saccEventsCutted.type(i) = "";
            saccEventsCutted.givenDuration(i) = NaN;
        end
    end
    eventDiff = length(newSaccadesTime) - length(newFixationsTime);
    %% Save trigger file
    writetable(fixEventsCutted,[triggerFolder,filesep,'rawTriggerFileFixations_Participant',participants(sub).name(end-7:end-4),'.csv']);
    save([EEGdataFolder,filesep,'preproc_',participants(sub).name(1:end-4),filesep,'alignmentMetadata',participants(sub).name(1:end-4),'.mat'],'differenceNumberMessages','residualsTimeline','linearDrift','intercept','eventDiff')

    writetable(saccEventsCutted,[triggerFolder,filesep,'rawTriggerFileSaccades_Participant',participants(sub).name(end-7:end-4),'.csv']);
end
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc