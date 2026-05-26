%% Cosine similarity computation
% This script was developed by Aitana Grasso-Cladera
%% History
% 20.04.26 - AGC = Developed the script for computing cosine similarities
% for all pairs or mass-univariate from Unfold. The scripts also averages
% the cosine similarity coefficients for subjects with more than one
% session, and do the grand average across all subjects. It implements
% Fishers' Z tranformation before average and saves the data for plotting.
%% What type of data?
power = 0;
ERP = 1;
%% Load the data
projectFolder = [];
dataFolder = [projectFolder,filesep,'Data'];
unfoldFolder = [dataFolder,filesep,'unfold'];
if power == 1
    load(unfoldFolder,filesep,'betasUnfoldALL_overlapFinalMatchingBlinksTF.mat');
elseif ERP == 1
    load(unfoldFolder,filesep,'betasUnfoldALL_overlapFinalMatchingBlinks.mat');
end
%% Compute cosine similarity between target mass-univariate for each 
% session and subject
pairs = {'BlinkOffset','Fixation';'BlinkOffset','BlinkOnset';'BlinkOffset',...
    'Saccade';'Fixation','BlinkOnset';'Fixation','Saccade';'BlinkOnset',...
    'Saccade'};

similarity = nan(size(blinkOn,1),size(blinkOn,2),size(pairs,1));
similarityTransformed = nan(size(blinkOn,1),size(blinkOn,2),size(pairs,1));

for j = 1:size(blinkOn,1)
    for ch = 1:size(blinkOn,2)
        % Get the data for every session and every channel
        blinkOffset = squeeze(mean(blinkOff(j,ch,:),1))'; %1
        fixation = squeeze(mean(fix(j,ch,:),1))'; %2
        blinkOnset = squeeze(mean(blinkOn(j,ch,:),1))'; %3
        saccade = squeeze(mean(sacc(j,ch,:),1))'; %4
        
        for i = 1:size(pairs,1)
            ERP1 = pairs{i,1};
            ERP2 = pairs{i,2};

            if contains(ERP1,'BlinkOffset')
                erpA = blinkOffset;
            elseif contains(ERP1,'Fixation')
                erpA = fixation;
            elseif contains(ERP1,'BlinkOnset')
                erpA = blinkOnset;
            elseif contains(ERP1,'Saccade')
                erpA = saccade;
            end

            if contains(ERP2,'BlinkOffset')
                erpB = blinkOffset;
            elseif contains(ERP2,'Fixation')
                erpB = fixation;
            elseif contains(ERP2,'BlinkOnset')
                erpB = blinkOnset;
            elseif contains(ERP2,'Saccade')
                erpB = saccade;
            end
            % Compute cosine similarity = (A·B) / (|A| * |B|)
            numerator   = dot(erpA, erpB);
            denominator = norm(erpA) * norm(erpB);

            if denominator == 0
                sim = NaN;
            else
                sim = numerator / denominator;
            end

            similarity(j,ch,i) = sim;

            similarityTransformed(j,ch,i) = 0.5 * (log(1+sim) - log(1-sim));
         end
    end
end
%% Average over sessions of the same subject
dataFolder = unfoldFolder;
%%
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.'|| contains(tmp(pId).name,'TF')
        continue
    else
        participants(inx).name = tmp(pId).name;
        participants(inx).folder = tmp(pId).folder;
        participants(inx).date = tmp(pId).date;
        inx = inx+1;
    end
end

clear tmp
%% Find which subjects have 2 sessions
% Get names from all files
subjectNumber = nan(1,length(participants)-1);
for i = 2:length(participants)
    subjectNumber(i) = str2double(participants(i).name(end-7:end-6));
end
subjectNumber(1) = [];
[uniqueVals, ~, groupIdx] = unique(subjectNumber,'stable');
nSubjects = numel(uniqueVals);
%% Take all the coefficients for each pair of comparison
blinkOff_fixation = similarityTransformed(:,:,1);
blinkOff_blinkOn = similarityTransformed(:,:,2);
blinkOff_saccade = similarityTransformed(:,:,3);
fixation_blinkOn = similarityTransformed(:,:,4);
fixation_saccade = similarityTransformed(:,:,5);
blinkOn_saccade = similarityTransformed(:,:,6);

final_blinkOff_fixation = nan(nSubjects,size(blinkOff,2));
final_blinkOff_blinkOn  = nan(nSubjects,size(blinkOn,2));
final_blinkOff_saccade = nan(nSubjects,size(fix,2));
final_fixation_blinkOn  = nan(nSubjects,size(sacc,2));
final_fixation_saccade  = nan(nSubjects,size(sacc,2));
final_blinkOn_saccade  = nan(nSubjects,size(sacc,2));

for s = 1:nSubjects
    idx = find(groupIdx == s);   % rows belonging to this subject
    
    if numel(idx) == 1
        % only one session = copy
        final_blinkOff_fixation(s,:) = blinkOff_fixation(idx,:);
        final_blinkOff_blinkOn(s,:)  = blinkOff_blinkOn(idx,:);
        final_blinkOff_saccade(s,:) = blinkOff_saccade(idx,:);
        final_fixation_blinkOn(s,:)  = fixation_blinkOn(idx,:);
        final_fixation_saccade(s,:)  = fixation_saccade(idx,:);
        final_blinkOn_saccade(s,:)  = blinkOn_saccade(idx,:);
    else
        % two (or more) sessions = average
        final_blinkOff_fixation(s,:) = mean(blinkOff_fixation(idx,:),1);
        final_blinkOff_blinkOn(s,:)  = mean(blinkOff_blinkOn(idx,:),1);
        final_blinkOff_saccade(s,:) = mean(blinkOff_saccade(idx,:),1);
        final_fixation_blinkOn(s,:)  = mean(fixation_blinkOn(idx,:),1);
        final_fixation_saccade(s,:)  = mean(fixation_saccade(idx,:),1);
        final_blinkOn_saccade(s,:)  = mean(blinkOn_saccade(idx,:),1);
    end
end
%% Average across all subjects
mean_blinkOff_fixation = nan(1,size(blinkOff,2));
mean_blinkOff_blinkOn  = nan(1,size(blinkOn,2));
mean_blinkOff_saccade = nan(1,size(fix,2));
mean_fixation_blinkOn  = nan(1,size(sacc,2));
mean_fixation_saccade  = nan(1,size(sacc,2));
mean_blinkOn_saccade  = nan(1,size(sacc,2));

for ch = 1:size(final_blinkOn_saccade,2)
    mean_blinkOff_fixation(ch) = mean(final_blinkOff_fixation(:,ch),'omitnan');
    mean_blinkOff_blinkOn(ch)  = mean(final_blinkOff_blinkOn(:,ch),'omitnan');
    mean_blinkOff_saccade(ch)  = mean(final_blinkOff_saccade(:,ch),'omitnan');
    mean_fixation_blinkOn(ch)  = mean(final_fixation_blinkOn(:,ch),'omitnan');
    mean_fixation_saccade(ch)  = mean(final_fixation_saccade(:,ch),'omitnan');
    mean_blinkOn_saccade(ch)   = mean(final_blinkOn_saccade(:,ch),'omitnan');
end
%% Transform back to cosine coefficients, this is the data we will use to 
% plot
ch = 11;
plot_blinkOff_fixation = (exp(2*mean_blinkOff_fixation(ch)) - 1) / ...
                         (exp(2*mean_blinkOff_fixation(ch)) + 1);
plot_blinkOff_blinkOn = (exp(2*mean_blinkOff_blinkOn(ch)) - 1) / ...
                        (exp(2*mean_blinkOff_blinkOn(ch)) + 1);
plot_blinkOff_saccade = (exp(2*mean_blinkOff_saccade(ch)) - 1) / ...
                        (exp(2*mean_blinkOff_saccade(ch)) + 1);
plot_fixation_blinkOn = (exp(2*mean_fixation_blinkOn(ch)) - 1) / ...
                        (exp(2*mean_fixation_blinkOn(ch)) + 1);
plot_fixation_saccade = (exp(2*mean_fixation_saccade(ch)) - 1) / ...
                        (exp(2*mean_fixation_saccade(ch)) + 1);
plot_blinkOn_saccade = (exp(2*mean_blinkOn_saccade(ch)) - 1) / ...
                       (exp(2*mean_blinkOn_saccade(ch)) + 1);
if power == 1
    save(fullfile(unfoldFolder,sprintf('cosineSimilarityPlotTF.mat')),...
        'plot_blinkOff_fixation','plot_blinkOff_blinkOn','plot_blinkOff_saccade',...
        'plot_fixation_blinkOn','plot_fixation_saccade','plot_blinkOn_saccade');
elseif ERP == 1
    save(fullfile(unfoldFolder,sprintf('cosineSimilarityPlotERP.mat')),...
        'plot_blinkOff_fixation','plot_blinkOff_blinkOn','plot_blinkOff_saccade',...
        'plot_fixation_blinkOn','plot_fixation_saccade','plot_blinkOn_saccade');
end