%% Unfold for Power on the Gamma Band
% This script was developed by Aitana Grasso-Cladera and Debora Nolte
%% History
% 19.04.26 - AGC = Developed the current script based on the ERP version.
%% Path definition
projectFolder = [];
dataFolder = [projectFolder,filesep,'Data'];
addpath('[]/eeglab2026.0.0');
workingFolder = [];
addpath(genpath('[]/unfold'))
init_unfold
addpath('[]/unfold/gramm')
%%
close
%%
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.'|| ~contains(tmp(pId).name,'preproc_')
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
cnt = 1;
for sub = 1:length(participants)
  
    participantFolder = [dataFolder,filesep,participants(sub).name];
    % add this new folder to the savedata path, so where the intermediate steps
    % will be saved
    savedata = [dataFolder,filesep,'unfold'];
    if ~exist(savedata,'dir')
        mkdir(savedata)
    end
    %% Load the data
    EEG = pop_loadset(sprintf('2a_cleanDataChannels_noRejectionAll_%s.set',participants(sub).name(9:end)),fullfile([dataFolder,filesep,participants(sub).name]));
    %% Interpolating missing channels
    % Get all channels that need to be interpolated
    EEGchan = pop_loadset(sprintf('1a_triggersFilteredAll_%s.set',participants(sub).name(9:end)),fullfile([dataFolder,filesep,participants(sub).name]));

    fullChanlocs = EEGchan.chanlocs; % Used for data cleaning and interpolation
    clear EEGchan
    EEG = pop_interp(EEG,fullChanlocs,'spherical');
    EEG = eeg_checkset(EEG);
    %% Prepare the event file
    % Remove boudary events (for now)
    boundaryEvents = false(size(EEG.event));
    for j = 1:length(EEG.event)
        if contains(EEG.event(j).type, 'boundary')
            boundaryEvents(j) = true;
        end
    end

    EEG.event(boundaryEvents) = [];

    % Remove invalid fixations and saccades
    invalidEvents = false(size(EEG.event));
    for i = 1:length(EEG.event)
        if contains(EEG.event(i).status, 'invalid')
            invalidEvents(i) = true;
        end
    end
    EEG.event(invalidEvents) = [];
    %% Create design matrix
    cfgDesign               = [];
    % all events we have are saccades, blinks and fixations:
    cfgDesign.eventtypes = {'blinkOnset','blinkOffset','saccadeOnset','fixationOnset'};
    cfgDesign.codingschema  = 'effects';
    cfgDesign.formula       = {'y ~ 1','y ~ 1','y ~ 1','y ~ 1'};
   % cfgDesign.channel       = 11;
    
    EEG                    = uf_designmat(EEG,cfgDesign);  
    %% Timeshift
    cfgTimeshift            = [];
    cfgTimeshift.timelimits = [-1 1];
    EEG                    = uf_timeexpandDesignmat(EEG,cfgTimeshift);
    %% Reject noisy segments
    % % we only do this, if we have noisy segments to reject
    if isfile(fullfile([dataFolder,filesep,participants(sub).name],sprintf('removed_intervalsCleanRawData_%s.mat',participants(sub).name(9:end))))
        load(fullfile([dataFolder,filesep,participants(sub).name],sprintf('removed_intervalsCleanRawData_%s.mat',participants(sub).name(9:end))));
        % Check if we are removing the blink completely (on + off), if not,
        % do it
        eventTypes = {EEG.event.type};
        onsetIdx  = find(strcmp(eventTypes,'blinkOnset'));
        offsetIdx = find(strcmp(eventTypes,'blinkOffset'));

        newWinrej = [];
        for i = 1:length(onsetIdx)
            onset_i   = onsetIdx(i);
            onset_lat = EEG.event(onset_i).latency;
            % Find first offset occurring after this onset
            offset_i = offsetIdx(find(offsetIdx > onset_i,1,'first'));
            if isempty(offset_i)
                continue
            end
            offset_lat = EEG.event(offset_i).latency;
            % Check whether onset or offset is already rejected
            onsetRejected  = any(onset_lat  >= tmprej(:,1) & onset_lat  <= tmprej(:,2));
            offsetRejected = any(offset_lat >= tmprej(:,1) & offset_lat <= tmprej(:,2));
            % If either one is rejected → reject the full blink interval
            if onsetRejected || offsetRejected
                newWinrej = [newWinrej; onset_lat offset_lat];
            end
        end
        % Merge with the matrix coming from clean raw data, sort them and
        % merge overlapping windows
        tmprej = [tmprej; newWinrej];
        if ~isempty(tmprej)
            tmprej = sortrows(tmprej);
            merged = tmprej(1,:);
            for i = 2:size(tmprej,1)
                if tmprej(i,1) <= merged(end,2)
                    merged(end,2) = max(merged(end,2), tmprej(i,2));
                else
                    merged = [merged; tmprej(i,:)];
                end
            end
            tmprej = merged;
        end

        EEG                = uf_continuousArtifactExclude(EEG,'winrej',tmprej);
    end
    %% Load detected noisy segments after clean raw data removal
    if isfile([participantFolder,filesep,sprintf('removed_IQRcontinuousFinalMatchingBlinks_%s.mat',participants(sub).name(9:end))])
       load([participantFolder,filesep,sprintf('removed_IQRcontinuousFinalMatchingBlinks_%s.mat',participants(sub).name(9:end))])
        EEG                = uf_continuousArtifactExclude(EEG,'winrej',bad_windowsTimes);
    end

    EEG1 = EEG;
   %% Compute TF (Morlet Wavelet based on Mike X Cohen)
   srate = 500;
   data  = EEG.data;        % channels x time
   [nChan,nData] = size(data);

   frex = logspace(log10(30),log10(80),20);
   nFreq = length(frex);
   nCycles = 7; % 3 is Peter's recommendation but looks very noisy
   wavetime = -1:1/srate:1;   % kernel support (just computational padding)
   nWave = length(wavetime);
   nConv = nWave + nData - 1;
   tf = zeros(nChan,nFreq,nData);
   tf_log = zeros(nChan,nFreq,nData);
   gamma_power = zeros(nChan,nData);

   for ch = 1:nChan
       dataX = fft(data(ch,:),nConv);
       for fi = 1:nFreq
           % Gaussian width determined by cycles
           s = nCycles/(2*pi*frex(fi));
           % Morlet wavelet
           wavelet = exp(2*1i*pi*frex(fi).*wavetime) ...
               .* exp(-wavetime.^2./(2*s^2));
           waveletX = fft(wavelet,nConv);
           waveletX = waveletX./max(waveletX);
           convRes = ifft(waveletX .* dataX);
           % Remove convolution edges
           convRes = convRes(floor(nWave/2)+1:end-floor(nWave/2));
           % COMPUTE POWER
           tf(ch,fi,:) = abs(convRes).^2;
       end
       % dB Normalization
       baseline = mean(tf(ch,:,:),3);
       tf_log(ch,:,:) = 10*log10(tf(ch,:,:)./baseline);

       % Average gamma band
       gamma_idx = frex >= 30 & frex <= 80;
       gamma_power(ch,:) = squeeze(mean(tf_log(ch,gamma_idx,:),2));
   end
   timeVector = (0:nData-1)/srate;
   %% Reeplace in structure
   EEG.data = gamma_power;
   EEG.times = timeVector;
   EEG = eeg_checkset(EEG);
   %% Fit the model
   % fit the model:
   % the parameters are taken from Gert et al., 2022
   cfgFit                  = [];
   cfgFit.precondition     = 1;
   cfgFit.lsmriterations   = 1500; % steps iterative solver should reach
   %cfgFit.channel          = 1; % all channels
    EEG                    = uf_glmfit(EEG, cfgFit);
    %% Make a massive uni-variate fit without de-convolution (Gert et al., 2022)
    EEGepoch = uf_epoch(EEG,'timelimits',cfgTimeshift.timelimits);

    EEGepoch = uf_glmfit_nodc(EEGepoch);
    %% Get the betas
    % results condensed in new structure
    ufresult                = uf_condense(EEG);
    ufresultEp              = uf_condense(EEGepoch);

    ufresultEp = uf_predictContinuous(ufresultEp); 
    ufresultEp = uf_addmarginal(ufresultEp);

    paramNames={ufresultEp.param.name};
    [~,paramBlinkOnset]=find(ismember(paramNames,'(Intercept)'));
    [~,paramBlinkOffset]=find(ismember(paramNames,'2_(Intercept)'));
    [~,paramSaccade]=find(ismember(paramNames,'3_(Intercept)'));
    [~,paramFixation]=find(ismember(paramNames,'4_(Intercept)'));

    sacc(cnt,:,:) = ufresultEp.beta(:,:,paramSaccade);
    blinkOn(cnt,:,:)  = ufresultEp.beta(:,:,paramBlinkOnset);
    blinkOff(cnt,:,:)  = ufresultEp.beta(:,:,paramBlinkOffset);
    fix(cnt,:,:)  = ufresultEp.beta(:,:,paramFixation);

    time = ufresultEp.times;

    cnt = cnt + 1; % to avoid having nans in the data

    %% Save the data
    save(fullfile(savedata,sprintf('betasUnfoldFinalMatchingBlinksTF_%s.mat',participants(sub).name)),'ufresult','ufresultEp','tf_log','-v7.3');
end

save(fullfile(saveData,sprintf('betasUnfoldALL_overlapFinalMatchingBlinksTF_.mat')),'sacc','time','fix','blinkOn','blinkOff');
%%
rmpath('[]/eeglab2025.0.0');
rmpath('[]/unfold')