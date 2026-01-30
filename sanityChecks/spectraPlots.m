%% this code plots power spectrum for all recordings

clear all;
close all;
clc;
%% 
basePath = '';
addpath('');
eeglab; close; %so the window doesnt pop up

folders = dir(fullfile(basePath, 'preproc_Participant*'));

targetChan = 'Cz';

% create a frequency axis (0–160 Hz, 0.5 Hz resolution) - this is needed bc
% all recordings have diff number of timepoints
common_hz = 0:0.5:160;

% storage
all_power_db = [];
participant_ids = {};

for i = 1:length(folders)
    if folders(i).isdir
        participantFolder = folders(i).name;
        folderPath = fullfile(basePath, participantFolder);

        % get participant ID
        tokens = regexp(participantFolder, 'Participant(\d+_\d+)', 'tokens');
  
        participantID = tokens{1}{1};

        % get file
        fileName = sprintf('5a_epochedFixationRejectionIQ_Participant%s.set', participantID);
        filePath = fullfile(folderPath, fileName);

        % Load EEG dataset
        EEG = pop_loadset('filename', fileName, 'filepath', folderPath);

        % get channel data
        chanData = EEG.data(chanIdx, :, :);

        % Compute FFT and Power
        fft_data = fft(chanData, [], 2) / EEG.pnts; % 2 - as the time is in the second dimention of the eeg data
        power = abs(fft_data).^2;
        power = power(:, 1:floor(EEG.pnts/2)+1, :); % keep only the positive frequencies
        mean_power = mean(power, 3); %Average across trials (trials are 3rd dimension of eeg data)
        mean_power_db = 10 * log10(mean_power);

        % Frequency axis
        hz = linspace(0, EEG.srate/2, size(mean_power_db, 2));

        % Interpolate to common frequency axis - bc all recordings have diff number of timepoints
        interp_power = interp1(hz, mean_power_db, common_hz, 'linear', 'extrap'); 

        % baseline correction
        % subtract the mean power across all frequencies
        baseline = mean(interp_power); 
        interp_power = interp_power - baseline;

        % Store
        all_power_db(i, :) = interp_power; 
        participant_ids{i} = participantID;

        fprintf('Processed %s (%s)\n', participantFolder, targetChan);
    end
end

% Plot all participants' Cz power spectra
figure('Name', sprintf('Power Spectrum at %s', targetChan)), clf
plot(common_hz, all_power_db', 'LineWidth', 1.2)
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title(sprintf('Power Spectrum at %s (All Participants)', targetChan))
axis([0 160 -50 50])       % fixed frequency & power scale
grid on
legend(participant_ids, 'Location', 'northeastoutside')
%% a check
disp(participant_ids')
%% to plot good and bad recordings separately
% for no IQ
%good = {'09_1','09_2','10_2','11_1','11_2','12_1','12_2','13_1','19_1','20_1','20_2','22_1','30_1','30_2','33_1','33_2','35_1','36_2'};

%for IQ
good =  {'09_1','09_2','10_2','11_1','11_2','12_1','12_2','13_1','16_1', '19_1','20_1','20_2','22_1','22_2', '30_1','30_2','33_1','33_2','35_1','36_2'};

participant_ids_str = string(participant_ids);

% find indices
goodIdx = ismember(participant_ids_str, good);
badIdx = ~goodIdx;

% Plot "good" participants 
figure('Name', 'Power Spectrum - "good"'), clf
plot(common_hz, all_power_db(goodIdx, :)', 'LineWidth', 1.2)
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('Cz Power Spectrum - "good"')
axis([0 160 -50 50])
grid on
legend(participant_ids_str(goodIdx), 'Location', 'northeastoutside')

% Plot "bad" participants
figure('Name', 'Power Spectrum - "bad"'), clf
plot(common_hz, all_power_db(badIdx, :)', 'LineWidth', 1.2)
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('Cz Power Spectrum - "bad"')
axis([0 160 -50 50])
grid on
legend(participant_ids_str(badIdx), 'Location', 'northeastoutside')


