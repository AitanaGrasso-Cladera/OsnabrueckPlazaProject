clear all;
close all;
clc;
%% loading the files, eeglab
basePath = '';
participantFolder = '';
fileName = ''; 

filePath = fullfile(basePath, participantFolder);

addpath(''); %eeglab
[ALLEEG, ~, CURRENTSET, ALLCOM] = eeglab;

EEG = pop_loadset('filename', fileName, 'filepath', filePath);

%% set-up
% choosing channel
chanLabels = {EEG.chanlocs.labels};
Cz = find(strcmpi(chanLabels,'Cz'));

times = EEG.times;       % time vector
win   = [140 200];      % window for N170/Vpp change for lpp

idx(1) = find(times >= win(1), 1, 'first'); % first sample >= 140 ms
idx(2) = find(times <= win(2), 1, 'last');  % last sample <= 200 ms

nTrials = size(EEG.data,3);
face_trials = [];
background_trials = []; 

for t = 1:nTrials % separate for faces and backgrounds
    evt = EEG.epoch(t).eventtype; 
    if any(strcmp(evt, 'face'))
        face_trials(end+1) = t;
    elseif any(strcmp(evt, 'background'))
        background_trials(end+1) = t;
    end
end

% baseline window (prestimulus)
%baselineWin = [-200 0];   % ms

%baseIdx(1) = find(times >= baselineWin(1), 1, 'first');
%baseIdx(2) = find(times <= baselineWin(2), 1, 'last');

%% local mean amplitude and latency

function [meanAmp, meanLat] = computeMean(EEG, chan, trialIdx, idx) %baseIdx)

    nTrials = length(trialIdx);
    meanAmp = zeros(nTrials,1);
    meanLat = zeros(nTrials,1);


    for i = 1:nTrials
        tr = trialIdx(i);

        % extract full epoch for baseline
        epoch = squeeze(EEG.data(chan, :, tr));

        % compute baseline offset
        %baseline = mean(epoch(baseIdx(1):baseIdx(2)));

        % baseline-corrected signal in analysis window
        sig_bc = epoch(idx(1):idx(2)); %- baseline;

        % mean amplitude
        meanAmp(i) = mean(sig_bc);

    end
end

[meanAmp_faces, meanLat_faces] = computeMean(EEG, Cz, face_trials, idx);
[meanAmp_bg, meanLat_bg]       = computeMean(EEG, Cz, background_trials, idx);

% means
meanAmp_faces_mean = mean(meanAmp_faces);
meanAmp_bg_mean = mean(meanAmp_bg);

% standard deviations
stdAmp_faces = std(meanAmp_faces);
stdAmp_bg = std(meanAmp_bg);

SummaryTable = table({'Faces'; 'Background'}, ...
    [meanAmp_faces_mean; meanAmp_bg_mean], ...
    [stdAmp_faces; stdAmp_bg], ...
    'VariableNames', {'Condition', ...
                      'MeanAmp','StdAmp'});
disp(SummaryTable);

%% local mean peak amplitude and latency
%{
% compute amplitude and latency
function [negAmp, negLat, posAmp, posLat] = computePeaks(EEG, chan, trialIdx, idx, baseIdx, times)

    nTrials = length(trialIdx);
    negAmp = zeros(nTrials,1);
    negLat = zeros(nTrials,1);
    posAmp = zeros(nTrials,1);
    posLat = zeros(nTrials,1);

    for i = 1:nTrials
        tr = trialIdx(i);

        % extract full epoch
        epoch = squeeze(EEG.data(chan, :, tr));

        % baseline correction
        baseline = mean(epoch(baseIdx(1):baseIdx(2)));
        sig_bc = epoch(idx(1):idx(2)) - baseline;

        % negative peak
        [negAmp(i), relIdx] = min(sig_bc);
        negLat(i) = times(idx(1) + relIdx - 1);

        % positive peak
        [posAmp(i), relIdx] = max(sig_bc);
        posLat(i) = times(idx(1) + relIdx - 1);
    end
end

[negAmp_faces, negLat_faces, posAmp_faces, posLat_faces] = ...
    computePeaks(EEG, Cz, face_trials, idx, baseIdx, times);

[negAmp_bg, negLat_bg, posAmp_bg, posLat_bg] = ...
    computePeaks(EEG, Cz, background_trials, idx, baseIdx, times);

% means
meanNegAmp_faces = mean(negAmp_faces);
meanNegLat_faces = mean(negLat_faces);
meanPosAmp_faces = mean(posAmp_faces);
meanPosLat_faces = mean(posLat_faces);

meanNegAmp_bg = mean(negAmp_bg);
meanNegLat_bg = mean(negLat_bg);
meanPosAmp_bg = mean(posAmp_bg);
meanPosLat_bg = mean(posLat_bg);

% standard deviations
stdNegAmp_faces = std(negAmp_faces);
stdNegLat_faces = std(negLat_faces);
stdPosAmp_faces = std(posAmp_faces);
stdPosLat_faces = std(posLat_faces);

stdNegAmp_bg = std(negAmp_bg);
stdNegLat_bg = std(negLat_bg);
stdPosAmp_bg = std(posAmp_bg);
stdPosLat_bg = std(posLat_bg);

% summary table with mean and std
SummaryTable = table({'Faces'; 'Background'}, ...
    [meanNegAmp_faces; meanNegAmp_bg], ...
    [stdNegAmp_faces; stdNegAmp_bg], ...
    [meanNegLat_faces; meanNegLat_bg], ...
    [stdNegLat_faces; stdNegLat_bg], ...
    [meanPosAmp_faces; meanPosAmp_bg], ...
    [stdPosAmp_faces; stdPosAmp_bg], ...
    [meanPosLat_faces; meanPosLat_bg], ...
    [stdPosLat_faces; stdPosLat_bg], ...
    'VariableNames', {'Condition', ...
                      'MeanNegAmp','StdNegAmp', ...
                      'MeanNegLat','StdNegLat', ...
                      'MeanPosAmp','StdPosAmp', ...
                      'MeanPosLat','StdPosLat'});
disp(SummaryTable);
%}

%% permutation t-test to see if the difference is significant (same time window)

% Helper function for permutation t-test
function pVal = permTtest(group1, group2, nPerm)
    if nargin < 3 % n of input arguments
        nPerm = 1000; 
    end
    n1 = length(group1);
    n2 = length(group2);
    
    allData = [group1; group2];
    labels  = [ones(n1,1); zeros(n2,1)];
    
    % t-test
    tObs = (mean(group1) - mean(group2)) / sqrt(var(group1)/n1 + var(group2)/n2);
    
    % Permutations
    tPerm = zeros(nPerm,1);
    for p = 1:nPerm
        permLabels = labels(randperm(length(labels)));
        g1 = allData(permLabels==1);
        g2 = allData(permLabels==0);
        tPerm(p) = (mean(g1) - mean(g2)) / sqrt(var(g1)/length(g1) + var(g2)/length(g2));
    end
    
    pVal = mean(abs(tPerm) >= abs(tObs));
end

%pValNeg = permTtest(negAmp_faces, negAmp_bg, 1000);
%pValPos = permTtest(posAmp_faces, posAmp_bg, 1000);

pValmean = permTtest(meanAmp_faces, meanAmp_bg, 1000);

pValmean

%% Plot ERPs 
face_data = EEG.data(Cz, :, face_trials);   % time × trials
%{
for t = 1:size(face_data, 3)
    baseline = mean(face_data(1, baseIdx(1):baseIdx(2), t));
    face_data(1, :, t) = face_data(1, :, t) - baseline;
end
%}
erp_face = mean(face_data, 3);

bg_data = EEG.data(Cz, :, background_trials);
%{
for t = 1:size(bg_data, 3)
    baseline = mean(bg_data(1, baseIdx(1):baseIdx(2), t));
    bg_data(1, :, t) = bg_data(1, :, t) - baseline;
end
%}

erp_bg = mean(bg_data, 3);

figure;
plot(times, erp_face, 'b', 'LineWidth', 2); hold on;
plot(times, erp_bg, 'r', 'LineWidth', 2);

xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title('Average ERP at Cz: Faces vs Background');
legend('Faces', 'Background');
grid on;

% line for fix/sacc onset
yl = ylim;
plot([0 0], yl, 'k--');
xlim([-300 500]);

%% ERPs with random subset of background trials (same number as face trials)

% Number of face trials
nFace = length(face_trials);

% random subset of background trials (same number as face trials)
rng('shuffle');  % optional: shuffle random seed
rand_bg_idx = randperm(length(background_trials), nFace); 
bg_trials_subset = background_trials(rand_bg_idx);

% Baseline-correct random subset of background trials
bg_data_random = EEG.data(Cz, :, bg_trials_subset);   % time × trials

%{
for t = 1:size(bg_data_random, 3)
    baseline = mean(bg_data_random(1, baseIdx(1):baseIdx(2), t));
    bg_data_random(1, :, t) = bg_data_random(1, :, t) - baseline;
end
%}
erp_bg_random = mean(bg_data_random, 3);

% Plot
figure;
plot(times, erp_face, 'b', 'LineWidth', 2); hold on;
plot(times, erp_bg_random, 'r', 'LineWidth', 2);

xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title('Average ERP at Cz: Faces vs Random Background Trials');
legend('Faces', 'Background (subset)');
grid on;

% Line for fixation onset
yl = ylim;
plot([0 0], yl, 'k--')

%% ERPimage plots

figure;
erpimage(face_data, [], times, 'Cz – Face trials', 10, 1, ...
         'yerplabel','\muV','erp','on','cbar','on');

figure;
erpimage(bg_data, [], times, 'Cz – Background trials', 10, 1, ...
         'yerplabel','\muV','erp','on','cbar','on');

%% Time-frequency analysis (ERSP)

% parameters: cycles = [3 0.8]; nfreq = 40, [10 40], 'timesout', 200

%baseline_window = [-200 0];

%
[ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef({face_data bg_data}, ...
    EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, [3 0.8], ...
    'alpha', 0.01, 'commonbase', 'on','freqs', [10 40], ... 
    'nfreqs', 40, 'baseline', [-200 0], 'basenorm', 'off', ...
    'timesout', 200, 'trialbase', 'off', ...
    'plotersp', 'off', 'plotitc', 'off');
%'alpha', 0.01 tests each time-frequency point at p<0.01, 
% third plot shows only statistically significant differences

figure;
% Face ERSP: 
subplot(1,3,1);
imagesc(times, freqs, real(ersp{1}));
axis xy; xlabel('Time (ms)'); ylabel('Frequency (Hz)');
title('Face ERSP'); clim([-1 1]); colorbar;
hold on; plot([0 0], ylim, 'k--', 'LineWidth', 2); hold off;

% Background ERSP: 
subplot(1,3,2);
imagesc(times, freqs, real(ersp{2}));
axis xy; xlabel('Time (ms)'); ylabel('Frequency (Hz)');
title('Background ERSP'); clim([-1 1]); colorbar;
hold on; plot([0 0], ylim, 'k--', 'LineWidth', 2); hold off;

% Difference ERSP: (also computed by the NEWTIMEF for p<=0,01)
subplot(1,3,3);
imagesc(times, freqs, real(ersp{3}));
axis xy; xlabel('Time (ms)'); ylabel('Frequency (Hz)');
title('Face − Background'); clim([-1 1]); colorbar;
hold on; plot([0 0], ylim, 'k--', 'LineWidth', 2); hold off;