%% Average of subjects with more than one session
% This script was developed by Aitana Grasso-Cladera and Debora Nolte
%% What data are you working on?
power = 1;
ERP = 0;
%% Set paths
projectFolder = [];
dataFolder = [projectFolder,filesep,'Data',filesep,'unfold'];
%%
tmp = dir(fullfile(dataFolder));
participants = [];
inx = 1;

for pId = 1:size(tmp,1)
    if tmp(pId).name(1) == '.'|| ~contains(tmp(pId).name,'NEW') || contains(tmp(pId).name,'TF')
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
%%
if power == 1
    load([dataFolder,filesep,'betasUnfoldALL_overlapFinalMatchingBlinksTF.mat'])
elseif ERP == 1
    load([dataFolder,filesep,'betasUnfoldALL_overlapFinalMatchingBlinks.mat'])
end

finalBlinkOff = nan(nSubjects,size(blinkOff,2),size(blinkOff,3));
finalBlinkOn  = nan(nSubjects,size(blinkOn,2),size(blinkOn,3));
finalFixation = nan(nSubjects,size(fix,2),size(fix,3));
finalSaccade  = nan(nSubjects,size(sacc,2),size(sacc,3));

for s = 1:nSubjects
    
    idx = find(groupIdx == s);   % rows belonging to this subject
    
    if numel(idx) == 1
        % only one session = copy
        finalBlinkOff(s,:,:) = blinkOff(idx,:,:);
        finalBlinkOn(s,:,:)  = blinkOn(idx,:,:);
        finalFixation(s,:,:) = fix(idx,:,:);
        finalSaccade(s,:,:)  = sacc(idx,:,:);
    else
        % two (or more) sessions = average
        finalBlinkOff(s,:,:) = mean(blinkOff(idx,:,:),1);
        finalBlinkOn(s,:,:)  = mean(blinkOn(idx,:,:),1);
        finalFixation(s,:,:) = mean(fix(idx,:,:),1);
        finalSaccade(s,:,:)  = mean(sacc(idx,:,:),1);
    end
end
if power == 1
    save([dataFolder,filesep,'averagedTF.mat'],'finalBlinkOff','finalBlinkOn', ...
        'finalFixation','finalSaccade')
elseif ERP == 1
    save([dataFolder,filesep,'averagedERP.mat'],'finalBlinkOff','finalBlinkOn', ...
        'finalFixation','finalSaccade')
end