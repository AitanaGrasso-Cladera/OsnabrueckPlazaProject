clear all; % clear workspace
clc; % clear command window
close all; % close plots etc.

%% Path definition
dataFolder = '';
triggerFolder = '';
addpath('');
eventFolder = '';
workingFolder = '';

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
for sub = 1:length(participants)
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

    % Filter the data parameters adapted from Czeszumski, 2023 (Hyperscanning
    % Maastricht)
    low_pass = 128;
    high_pass = 10;
    EEG = pop_eegfiltnew(EEG, high_pass, []); % 10 is the lower edge
    EEG = pop_eegfiltnew(EEG, [], low_pass); % 128 is the upper edge

    %% Add central and right channels
    channels_to_sum = {'Cz', 'C3', 'C4', 'R1', 'R2', 'R3', 'R4'}; % Define the channels to be summed
    idx = find(ismember({EEG.chanlocs.labels}, channels_to_sum)); % Find their ids
    summed_signal = sum(EEG.data(idx, :, :), 1); % Sum the values

    EEG.data(end+1, :, :) = summed_signal;

    new_label = 'summed_channel'; % Define the name of the channel
    EEG.chanlocs(end+1).labels = 'summed_channel';

    EEG.nbchan = size(EEG.data, 1);

    %% Import the events file, find external triggers
    % Load the events file
    events = readtable([eventFolder, filesep,'events_',participants(sub).name(end-7:end-4)], 'TextType', 'string');
    events = renamevars(events, 'timestamp_ns_', 'timestamp_ns'); % Change the name of the variable
    events.timestamp_s = events.timestamp_ns / 1e9; % Turn ns to s
    events([1, 2, end], :) = []; % Delete the first two and the last rows
    events.timestamp_s = events.timestamp_s - min(events.timestamp_s); % Make the starting point zero
    trigger_types = ["A", "B", "C", "D"]; % Define the external trigger event names
    trigger_events = events(ismember(events.name, trigger_types), :); % Select the external triggers in the events file

    % Transform timestamps into sampling points
    latencies = trigger_events.timestamp_s;
    latencies = latencies * EEG.srate;

    % Define the search window
    window_ms = 300;
    window_pts = round((window_ms / 1000) * EEG.srate); % turn it into sample points

    % Define the channel
    channel_name = 'summed_channel';

    % Find the index of the channel of interest
    chan_idx = find(strcmp({EEG.chanlocs.labels}, channel_name));
    if isempty(chan_idx)
        error('summed_channel not found in EEG data.');
    end

    % Detect the events 
    for i = 1:length(latencies)
    
        start_latency = latencies(i) - window_pts; % Define start of the search window
        end_latency = latencies(i) + window_pts; % Define end of the search window

        segment = EEG.data(chan_idx, start_latency:end_latency); % Define the search window in the data

        segment = segment - mean(segment); % Baseline correction: subtract mean of the segment

        [peaks, locs] = findpeaks(segment); % Find positive peaks

        % If at least two peaks are found
        if length(peaks) >= 2
            % Sort peaks by their amplitude
            [sorted_peaks, sort_idx] = sort(peaks, 'descend'); % Sort in descending order (amplitude)
            top_locs = locs(sort_idx(1:2)); % Get the location of the two highest peaks
                
            % Sort the locations of the top peaks by time order (ascending)
            [sorted_locs, time_sort_idx] = sort(top_locs); % Sort the peak locations by time
            top_locs = sorted_locs; % Re-assign the locations in time order

            if top_locs(1) > top_locs(2)
                top_loc = top_loc(2);
            else
                top_loc = top_locs(1);
            end

        elseif length(peaks) == 1
           top_locs = locs;
           top_loc = locs; % If only one peak is found, take it
        else
           top_locs = []; % No peaks found
           top_loc = [];
        end



        if ~isempty(top_loc)
            % Add the top peak as a new event
        
            peak_latency = start_latency + top_loc - 1; % Convert to EEG time
                
            % Create new event with unique label (based on time order)
            EEG.event(end+1).type = trigger_events.name(i) + '_on';
            EEG.event(end).latency = peak_latency;
            EEG = eeg_checkset(EEG, 'eventconsistency');
            EEG = eeg_checkset(EEG, 'makeur'); % make urevent field
        end
    end
    EEG = pop_saveset(EEG, 'filename', [participants(sub).name(1:end-4), '_with_events'],'filepath',fullfile(participantFolder));
end

