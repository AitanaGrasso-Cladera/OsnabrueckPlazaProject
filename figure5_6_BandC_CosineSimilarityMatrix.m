%% Figure 5 and 6 B and C. Time-series (5) and Power (6)
%% Paths
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%% Which data are you working on?
ERP = 1;
POWER = 0;
%% Load the data
dataFolder = [projectFolder,filesep,'Data',filesep,'unfold'];
if ERP == 1
    load([dataFolder,filesep,'betasUnfoldALL_overlapFinalAllEvents.mat'])
elseif POWER == 1
    load([dataFolder,filesep,'betasUnfoldALL_overlapFinalAllEventsTF.mat'])
end
%% Compute cosine similarity between target ERPs for each session and
% subject
pairs = {'BlinkOffset','Fixation';'BlinkOffset','BlinkOnset';'BlinkOffset',...
    'Saccade';'Fixation','BlinkOnset';'Fixation','Saccade';'BlinkOnset','Saccade'};

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
        % only one session → copy
        final_blinkOff_fixation(s,:) = blinkOff_fixation(idx,:);
        final_blinkOff_blinkOn(s,:)  = blinkOff_blinkOn(idx,:);
        final_blinkOff_saccade(s,:) = blinkOff_saccade(idx,:);
        final_fixation_blinkOn(s,:)  = fixation_blinkOn(idx,:);
        final_fixation_saccade(s,:)  = fixation_saccade(idx,:);
        final_blinkOn_saccade(s,:)  = blinkOn_saccade(idx,:);
    else
        % two (or more) sessions → average
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
%% Transform back to cosine coefficients
ch = 11;
plot_blinkOff_fixation = (exp(2*mean_blinkOff_fixation) - 1) / ...
                         (exp(2*mean_blinkOff_fixation) + 1);

plot_blinkOff_blinkOn = (exp(2*mean_blinkOff_blinkOn) - 1) / ...
                        (exp(2*mean_blinkOff_blinkOn) + 1);

plot_blinkOff_saccade = (exp(2*mean_blinkOff_saccade) - 1) / ...
                        (exp(2*mean_blinkOff_saccade) + 1);

plot_fixation_blinkOn = (exp(2*mean_fixation_blinkOn) - 1) / ...
                        (exp(2*mean_fixation_blinkOn) + 1);

plot_fixation_saccade = (exp(2*mean_fixation_saccade) - 1) / ...
                        (exp(2*mean_fixation_saccade) + 1);

plot_blinkOn_saccade = (exp(2*mean_blinkOn_saccade) - 1) / ...
                       (exp(2*mean_blinkOn_saccade) + 1);
%% plot
pairNames = {'blinkOffset_vs_fixation'   plot_blinkOff_fixation(ch)
    'blinkOffset_vs_blinkOnset' plot_blinkOff_blinkOn(ch)   
    'blinkOffset_vs_saccade'    plot_blinkOff_saccade(ch)
    'fixation_vs_blinkOnset'    plot_fixation_blinkOn(ch)
    'fixation_vs_saccade'       plot_fixation_saccade(ch)
    'blinkOnset_vs_saccade'     plot_blinkOn_saccade(ch)};

pairStr = pairNames(:,1);
pairVal = cell2mat(pairNames(:,2));

% Manual order to highlight the significant pairs
allLabels = cellfun(@(s) regexp(s, '(.*?)_vs_(.*)', 'tokens'), ...
                    pairStr, 'UniformOutput', false);

labels = {'blinkOffset','fixation','blinkOnset','saccade'};
n = numel(labels);
S = eye(n);
for k = 1:numel(pairVal)
    a = allLabels{k,1}{1,1}{1,1};
    b = allLabels{k,1}{1,1}{1,2};
    i = find(strcmp(labels,a));
    j = find(strcmp(labels,b));
    S(i,j) = pairVal(k);
    S(j,i) = pairVal(k);
end

displayLabels = {'Blink Offset','Fixation Onset','Blink Onset','Saccade Onset'};
cmapERP = magma(70);
fontSize = 16;

ax1 = figure;
imagesc(S)
hold on
axis equal tight
set(gca,'XTick',1:n,'XTickLabel',displayLabels,'XTickLabelRotation',45,...
        'YTick',1:n,'YTickLabel',displayLabels,'YDir','normal');
fontsize(fontSize,'points')
baseMap = magma(256);
idx = round(linspace(80,220,256));
myCmap = baseMap(idx,:);
colormap(myCmap)
clim([-1 1])
c = colorbar;
c.Label.String = 'Cosine similarity';
c.Ticks = [-1 0 1];
c.FontSize = fontSize;

textStrings = num2str(S(:),'%0.2f');
textStrings = strtrim(cellstr(textStrings));

[x,y] = meshgrid(1:n);

hText = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center',...
    'FontSize',fontSize,'FontWeight','bold');

midpoint = 0;

for k = 1:numel(S)
    if S(k) > midpoint
        set(hText(k),'Color','k')
    else
        set(hText(k),'Color','w')
    end
end

hold on
xline(2.5,'k','LineWidth',2)
yline(2.5,'k','LineWidth',2)
set(gca,'fontsize',fontSize)
pos1 = ax1.Position;
%%

if ERP == 1
    exportgraphics(gcf, fullfile(cd,'Figure 5B.png'),'Resolution', 600);
    save([dataFolder,filesep,'finalDatasetCosineERP.mat'],'plot_blinkOff_blinkOn','plot_blinkOff_saccade',...
        'plot_blinkOff_fixation','plot_fixation_blinkOn','plot_fixation_saccade','plot_blinkOn_saccade')
elseif POWER == 1
    exportgraphics(gcf, fullfile(cd,'Figure 6B.png'),'Resolution', 600);
    save([dataFolder,filesep,'finalDatasetCosineTF.mat'],'plot_blinkOff_blinkOn','plot_blinkOff_saccade',...
        'plot_blinkOff_fixation','plot_fixation_blinkOn','plot_fixation_saccade','plot_blinkOn_saccade')
end
close



