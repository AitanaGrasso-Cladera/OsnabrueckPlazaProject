%% Figure 4. Saccade- and blink-offset ERPs at Cz, sorted by saccade and
% blink duration, respectively, for individual participants.
%%
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%% Load data sorted data
load('binnedDataSacc.mat')
load('binnedDataBlink.mat')
%% Compute AIC difference
allEvents = load('betasUnfoldALL_overlapFinalAllEvents.mat');
offsetEvents = load('betasUnfoldALL_overlapFinalEndEvents.mat');
onsetEvents = load('betasUnfoldALL_overlapFinalStartEvents.mat');
middleEvents = load('betasUnfoldALL_overlapFinalMiddleEvents.mat');
% Find subjects with more than one session and average them
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
% Find which subjects have 2 sessions 
subjectNumber = nan(1,length(participants)-1);
for i = 2:length(participants)
    subjectNumber(i) = str2double(participants(i).name(end-7:end-6));
end
subjectNumber(1) = [];
[uniqueVals, ~, groupIdx] = unique(subjectNumber,'stable');
nSubjects = numel(uniqueVals);
for s = 1:nSubjects
    idx = find(groupIdx == s);   % rows belonging to this subject
    if numel(idx) == 1
        % only one session → copy
        AICallEvents(s) = allEvents.allAIC(s);
        AICoffsetEvents(s)  = offsetEvents.allAIC(s);
        AIConsetEvents(s)  = onsetEvents.allAIC(s);
        AICmiddleEvents(s)  = middleEvents.allAIC(s);
    else
        % two (or more) sessions → average
        AICallEvents(s) = mean(allEvents.allAIC(idx));
        AICoffsetEvents(s)  = mean(offsetEvents.allAIC(idx));
        AIConsetEvents(s) = mean(onsetEvents.allAIC(idx));
        AICmiddleEvents(s)  = mean(middleEvents.allAIC(idx));
    end
end
allAICS = [AICmiddleEvents',AIConsetEvents',AICoffsetEvents',AICallEvents'];
averageAIC = mean(allAICS,2);
AICdiffMatrix = [AICmiddleEvents'-averageAIC,AIConsetEvents'-averageAIC,AICoffsetEvents'-averageAIC,AICallEvents'-averageAIC];
%%
subjects = [2,5,7,13];
titles = {'Subject 1 - Session 2','Subject 3 - Session 2',...
    'Subject 4 - Session 2','Subject 12 - Session 2'};
figure(figure('Renderer', 'painters', 'Position', [10 10 1000 800]))
t = tiledlayout(5,2,'TileSpacing', 'Compact', 'Padding', 'Compact');
sacc = 1:2:8;
for i = 1:length(subjects)
    if i == 1
        ax1 = nexttile(sacc(i));
    end
    ax = nexttile(sacc(i));
    f = imagesc(times, 1:size(binnedDataSaccAll{1,subjects(i)}{11}, 2), ...
        binnedDataSaccAll{1,subjects(i)}{11}');
    colormap(ax,magma(5000))
    xline(0, 'k-', 'LineWidth', 2);
    col = colorbar;
    col.Label.String = 'μV';
    maxVal = max(abs(binnedDataSaccAll{1,subjects(i)}{11}(:)));
    clim([-maxVal maxVal])
    xlabel('Time (ms)','FontSize',14);
    ylabel('Bins of Trials','FontSize',14);
    hold on;
    nBins = size(binnedDataSaccAll{1,subjects(i)}{11}, 2);
    scatter(saccadeOnsetsALL{1,subjects(i)}, 1:nBins, 14, 'w', 'filled');
    hold off;
    title(titles{i}, 'Interpreter', 'none');
    set(gca,'fontsize',14)
end

blink = 2:2:8;
for i = 1:length(subjects)
    if i == 1
        ax2 = nexttile(blink(i));
    end
    ax = nexttile(blink(i));
    f = imagesc(times, 1:size(binnedDataBlinkALL{1,subjects(i)}, 2), ...
        binnedDataBlinkALL{1,subjects(i)}');
    colormap(ax,magma(5000))
    xline(0, 'k-', 'LineWidth', 2);
    col = colorbar;
    col.Label.String = 'μV';
    maxVal = max(abs(binnedDataBlinkALL{1,subjects(i)}(:)));
    clim([-maxVal maxVal])
    xlabel('Time (ms)','FontSize',14);
    ylabel('Bins of Trials','FontSize',14);
    hold on;
    nBins = size(binnedDataBlinkALL{1,subjects(i)}, 2);
    scatter(blinkOnsetsALL{1,subjects(i)}, 1:nBins, 14, 'w', 'filled');
    hold off;
    title(titles{i}, 'Interpreter', 'none');
    set(gca,'fontsize',14)
end

ax3 = nexttile([1 2]);
imagesc(AICdiffMatrix')
colormap(magma(20))
c = colorbar;
c.Label.String = {'AIC differece'};
c.Label.FontSize = 12;

xlabel('Subject')
ylabel('Model Events')
yticks(1:4)
yticklabels({'Onset','Middle','Offset','Onset+Offset'})

set(gca,'FontSize',12,'TickDir','out','Box','off')
axis tight
grid on

pos1 = ax1.Position;
pos2 = ax2.Position;
pos3 = ax3.Position;
fontSize = 16;
annotation('textbox', [pos1(1)-0.05, pos1(2)+pos1(4), 0.03, 0.03], ...
    'String','A','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')
annotation('textbox', [pos2(1)-0.05, pos2(2)+pos2(4), 0.03, 0.03], ...
    'String','B','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')
annotation('textbox', [pos3(1)-0.05, pos3(2)+pos3(4), 0.03, 0.03], ...
    'String','C','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

exportgraphics(gcf, fullfile(figureFolder,'Figure4.png'), ...
    'Resolution', 600);
close