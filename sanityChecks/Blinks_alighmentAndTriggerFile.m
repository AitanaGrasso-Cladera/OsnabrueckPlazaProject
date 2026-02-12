%% A copy of "alighmentAndTriggerFile.m" and a part of  "OsnaPlaza_preprocessing.m" scripts adjusted for blinks (onsets + offsets) used to create trigger files %%

clear all; %workspace
clc; % comman window
close all; %plots

%% getting files
%
% load eeglab and xdf plugin
addpath('');
addpath('');
eeglab; close; %so the window doesnt pop up

dataFolder = '';

% List participant folders
tmp = dir(fullfile(dataFolder, 'P*'));

participants = struct();
inx = 1;

for k = 1:numel(tmp)

    % Skip non-folders and hidden folders
    if ~tmp(k).isdir || tmp(k).name(1) == '.'
        continue
    end

    participants(inx).folder = fullfile(tmp(k).folder, tmp(k).name);
    participants(inx).name     = tmp(k).name;   
    participants(inx).id     = tmp(k).name(2:end);
    inx = inx + 1;
end

% Load blinks, events, xdf
for sub = 1:numel(participants)
    participantFolder = participants(sub).folder;
    % XDF file
    xdfFile = dir(fullfile(participantFolder, '*.xdf'));
    if isempty(xdfFile)
        error('No XDF file found in %s', participantFolder);
    end
    xdfPath = fullfile(xdfFile(1).folder, xdfFile(1).name);
    streams = load_xdf(xdfPath);
    % Find and save the RESTApi streams
    for i = 1:length(streams)
        if contains(streams{1,i}.info.name,'RESTApi')
            restAPIdata = streams{1,i};
        end
    end
    % Load events
    eventsFile = fullfile(participantFolder, 'events.csv');
    if ~exist(eventsFile, 'file')
        error('Missing events.csv in %s', participantFolder);
    end
    eventsPL = readtable(eventsFile, 'TextType', 'string');
    
    % Load blinks
    blinksFile = fullfile(participantFolder, 'blinks.csv');
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
        %% Get the ET events between start and end of EEG recording
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
        %% Get the ET events between start and end of EEG recording
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
    %% Create the trigger file 
    % read blinks
    blinkEvents = readtable(blinksFile, 'TextType', 'string');
    % Cut the blink data between EEG start and end
    lowerLimitETblink = min(find(blinkEvents.startTimestamp_ns_ >= eventsPL.timestamp_ns_(EEGstartMessage)));
    upperLimitETblink = max(find(blinkEvents.startTimestamp_ns_ <= eventsPL.timestamp_ns_(EEGendMessage)));

    blinkEvents = blinkEvents(lowerLimitETblink:upperLimitETblink,:);
    % Set ET timestamps from blinks to seconds
    blinkEventsStart = blinkEvents.startTimestamp_ns_/1e9;
    blinkEventsEnd = blinkEvents.endTimestamp_ns_/1e9; %for end ts
    % Set the first timestamp to zero
    blinkEventsTime = blinkEventsStart - firstETtimes(1);
    blinkEventsTimeEnd = blinkEventsEnd - firstETtimes(1); %for end ts
    % Round for precision
    blinkTimes = round(blinkEventsTime,3);
    blinkTimesEnd = round(blinkEventsTimeEnd,3); %for end ts
    %% Create an interporalted timeline
    fs = 1000;
    interpTime = linspace(EEGtimeInSeconds(1),EEGtimeInSeconds(end),(EEGtimeInSeconds(end)-EEGtimeInSeconds(1))*fs+1);
    interpTime = round(interpTime,3); % for precision
    %% Find blinks times on the inteprolated time
    actualBlinks = nan(size(blinkTimes));
    for i = 1:length(blinkTimes)
        if blinkTimes(i) <= interpTime(end)
            actualBlinks(i) = find(interpTime == blinkTimes(i));

        else
            actualBlinks(i) = NaN;
        end
    end
    actualBlinks(isnan(actualBlinks)) = [];
    %% Find blinks END times on the inteprolated time
    actualBlinksEnd = nan(size(blinkTimesEnd));
    for i = 1:length(blinkTimesEnd)
        if blinkTimesEnd(i) <= interpTime(end)
            actualBlinksEnd(i) = find(interpTime == blinkTimesEnd(i));

        else
            actualBlinksEnd(i) = NaN;
        end
    end
    actualBlinksEnd(isnan(actualBlinksEnd)) = [];
    %% Add the drift factor to every sample
    driftValues = intercept + linspace(0,linearDrift,length(interpTime));
    timeline = interpTime + driftValues;
    %% Select timepoints that correspond to blinks
    newBlinksTime = timeline(actualBlinks);
    newBlinksTimeEnd = timeline(actualBlinksEnd); % for end ts
    % Now we add the first EEG timestamps, to get the data in EEG latencies
    finalBlinksTime = (newBlinksTime*500)+1;
    finalBlinksTimeEnd = (newBlinksTimeEnd*500)+1; % for end ts
    %% Organize trigger files
    nRows = length(finalBlinksTimeEnd);  % number of rows % use finalBlinkTimeEnd size due to it being shorter after being cutted for EEG events (so there are no rows with onset and no offset time)
    % Preallocate each column
    blinkEventsCutted = table( ...
    nan(nRows,1), ... % blinkStart
    nan(nRows,1), ... % latency_onset 
    nan(nRows,1), ... % blinkEnd
    nan(nRows,1), ... % latency_offset 
    nan(nRows,1), ... % type
    nan(nRows,1), ... % givenDuration    
    'VariableNames', { ...
        'blinkStart', ...
        'latencyOnset', ...
        'blinkEnd', ...
        'latencyOffset', ...
        'type', ...
        'givenDuration'});

    for i = 1:length(finalBlinksTimeEnd) %closestIndex % use finalBlinkTimeEnd
            % Asign and save the relevant data
            blinkEventsCutted.blinkStart(i) = newBlinksTime(i);
            blinkEventsCutted.latencyOnset(i) = finalBlinksTime(i);
            blinkEventsCutted.blinkEnd(i) = newBlinksTimeEnd(i);
            blinkEventsCutted.latencyOffset(i) = finalBlinksTimeEnd(i);
            blinkEventsCutted.type(i) = NaN;
            blinkEventsCutted.givenDuration(i) = blinkEvents.duration_ms_(i);
    end
    %% Save trigger file
    writetable(blinkEventsCutted,[participantFolder,filesep,'rawTriggerFileBlinks_Participant',participants(sub).id,'.csv']);
    outputFile = fullfile(participantFolder, ['rawTriggerFileBlinks_Participant', participants(sub).id, '.csv']);
    writetable(blinkEventsCutted, outputFile);
    metadataFile = fullfile(participantFolder, ['alignmentMetadata_', participants(sub).id, '.mat']);
    %save(metadataFile, 'differenceNumberMessages', 'residualsTimeline', 'linearDrift', 'intercept');
  
    %save([participantFolder,filesep,'P',participants(sub).name(1:end-4),filesep,'alignmentMetadata',participants(sub).name(1:end-4),'.mat'],'differenceNumberMessages','residualsTimeline','linearDrift','intercept')
    

    %% from  "OsnaPlaza_preprocessing.m" script to create a clean trigger file
    hardwareTrigger = readtable([participantFolder,filesep,'hardwareTriggers_', participants(sub).id,'.csv'], 'TextType', 'string');
    % Arrange the trigger file
    data = readtable([participantFolder,filesep,'rawTriggerFileBlinks_Participant',participants(sub).id,'.csv'], 'TextType', 'string'); % Assuming the file has headers
        % Remove if there is any nan event
    nonValidTrial = find(isnan(data.givenDuration));
    if ~isempty(nonValidTrial)
            data(nonValidTrial,:) = [];
    end
        % Incorporate hardware triggers into the blinks trigger file
    newTabCol = nan(height(hardwareTrigger), 1);
    hardwareTrigger.("givenDuration") = newTabCol;
    hardwareTrigger.("blinkStart") = newTabCol;
    hardwareTrigger.("blinkEnd") = newTabCol;
    hardwareTrigger.("latencyOffset") = newTabCol;

    % rename latency column in hardware triggers to match blinks trigger
    % file and sort the data
    hardwareTrigger = renamevars(hardwareTrigger, "latency", "latencyOnset");

    data = [data;hardwareTrigger];
    data = sortrows(data,"latencyOnset");
        % Save the clean trigger file
    writetable(data,[participantFolder,filesep,'clean_','Participant',participants(sub).id,'_blinkTriggers.csv']);
  
end
clear all; % clear workspace
clc; % clear command window
close all; % close plots etc
