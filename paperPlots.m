%% Figure 3: Distribution of blink and saccade durations
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
triggerFolder = [projectFolder,filesep,'triggerFiles'];
addpath(genpath('[]/eeglab2026.0.0'));
workingFolder = [];
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];

tmp = dir(fullfile(rawDataFolder));
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
% For every session, extract the duration of each valid saccade and blink
for sub = 1:length(participants)
    participantFolder = [rawDataFolder,filesep,'preproc_',participants(sub).name(1:end-4)];

    if ~exist(participantFolder,'dir')
        mkdir(participantFolder)
    end
    dataSaccades = readtable([triggerFolder,filesep,'rawTriggerFileSaccadesNew_',participants(sub).name(1:end-4),'.csv'], 'TextType', 'string'); % Assuming the file has headers
    dataBlinks = readtable([triggerFolder,filesep,'blinkOnOff',filesep,'clean_',participants(sub).name(1:end-4),'_blinkTriggers.csv'], 'TextType', 'string'); % Assuming the file has headers

    saccStatus = dataSaccades.type;
    saccDuration = dataSaccades.givenDuration;
    for i = 1:size(saccStatus,1)
        if ismissing(saccStatus(i)) ||  contains(saccStatus(i),'invalid')
            continue
        end
        saccDur{sub}(i) = saccDuration(i);
    end
    blinkDuration = dataBlinks.givenDuration;
    for i = 1:size(blinkDuration,1)
        if isnan(blinkDuration(i))
            continue
        end
        blinkDur{sub}(i) = blinkDuration(i);
    end
end
% Plot into histograms
colors = magma(50);
binEdges = linspace(0, 550, 100);
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

nSub = length(blinkDur);

allCounts = zeros(nSub, length(binCenters));

figure;
% First, blink duration
ax1 = subplot(2,1,1);
hold on;
% Subject histograms (bars)
for sub = 1:nSub
    data = blinkDur{sub};
    data(data == 0) = [];
    % histogram per subject (normalized)
    counts = histcounts(data, binEdges, 'Normalization', 'probability');
    counts = smoothdata(counts, 'gaussian', 8);
    allCounts(sub,:) = counts;
    bar(binCenters, counts, ...
        'FaceAlpha', 0.12, ...
        'EdgeColor', 'none', ...
        'FaceColor', colors(25,:));
end
% Mean histogram
meanCounts = mean(allCounts, 1);
plot(binCenters, meanCounts, ...
    'k', ...
    'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Frequency');
grid on;
hold off;
a = gca;
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
set(gca,'fontsize',14)
% Now, saccade duration
ax2 = subplot(2,1,2);
binEdges = linspace(0, 500, 100);
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

nSub = length(blinkDur);

allCounts = zeros(nSub, length(binCenters));

hold on
for sub = 1:nSub
    data = saccDur{sub};
    data(data == 0) = [];
    % histogram per subject (normalized)
    counts = histcounts(data, binEdges, 'Normalization', 'probability');
    counts = smoothdata(counts, 'gaussian', 8);
    allCounts(sub,:) = counts;
    bar(binCenters, counts, ...
        'FaceAlpha', 0.12, ...
        'EdgeColor', 'none', ...
        'FaceColor', colors(38,:));
end

% Mean histogram
meanCounts = mean(allCounts, 1);
plot(binCenters, meanCounts, ...
    'k', ...
    'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Frequency');
grid on;
hold off;
a = gca;
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
set(gca,'fontsize',14)

pos1 = ax1.Position;
pos2 = ax2.Position;
fontSize = 16;
annotation('textbox', [pos1(1)-0.07, pos1(2)+pos1(4), 0.03, 0.03], ...
    'String','A','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')
annotation('textbox', [pos2(1)-0.07, pos2(2)+pos2(4), 0.03, 0.03], ...
    'String','B','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

exportgraphics(gcf, fullfile(figureFolder,'F3eventDurationDistribution.png'),'Resolution', 600);
close
%% Figure 4: Mass-univariate results for ERPs, and cross correlation
% Load data from ERPs after averaging sessions from same subject
dataFolder = [projectFolder,filesep,'Data',filesep,'unfold'];
load(fullfile(dataFolder,filesep,'averagedERP.mat'))
% Load cross correlation results
load(fullfile(dataFolder,filesep,'cosineSimilarityPlotERP.mat'))
% Prepare
pairNames = {'blinkOffset_vs_fixation'   plot_blinkOff_fixation
    'blinkOffset_vs_blinkOnset' plot_blinkOff_blinkOn   
    'blinkOffset_vs_saccade'    plot_blinkOff_saccade
    'fixation_vs_blinkOnset'    plot_fixation_blinkOn
    'fixation_vs_saccade'       plot_fixation_saccade
    'blinkOnset_vs_saccade'     plot_blinkOn_saccade};
pairStr = pairNames(:,1);
pairVal = cell2mat(pairNames(:,2));
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
% Compute confidence interval (95%)
ciBlinkOff = compute_CI(finalBlinkOff,.05);
ciBlinkOn = compute_CI(finalBlinkOn,.05);
ciSacc = compute_CI(finalSaccade,.05);
ciFix = compute_CI(finalFixation,.05);
fontSize = 22;

figure(figure('Renderer', 'painters', 'Position', [10 10 1900 1200]))
% First, we plot the mass univariated results
t = tiledlayout(1,2);
t.TileSpacing = "compact";
t.Padding = "compact";
ax1 = nexttile;
hold (ax1,'on')
hold on
[~,~] = plotCI(squeeze(mean(finalSaccade(:,11,:),1)),squeeze(ciSacc(11,:,:)),time,cmapERP(35,:),fontSize,'-');
[~,~] = plotCI(squeeze(mean(finalBlinkOn(:,11,:),1)),squeeze(ciBlinkOn(11,:,:)),time,cmapERP(35,:),fontSize,'--');
[~,~] = plotCI(squeeze(mean(finalFixation(:,11,:),1)),squeeze(ciFix(11,:,:)),time,cmapERP(55,:),fontSize,'-');
[~,~] = plotCI(squeeze(mean(finalBlinkOff(:,11,:),1)),squeeze(ciBlinkOff(11,:,:)),time,cmapERP(55,:),fontSize,'--');
xlim([min(time) max(time)])
legend('Saccade Onset','','Blink Onset','','Fixation Onset','','Blink Offset')
v = vline(0,'-');
v.Color = [0 0 0 0.8];
v.LineWidth = 1.5;
h = hline(0,'-');
h.Color = [0 0 0 0.8];
h.LineWidth = 1.5;
xlabel('Time (s)')
ylabel('Amplitude (µV)')
set(gca,'fontsize',fontSize)
% Now we plot the results from the cross correlation
ax2 = nexttile;
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
pos2 = ax2.Position;

annotation('textbox', [pos1(1)-0.03, pos1(2)+pos1(4), 0.03, 0.03], ...
    'String','A','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

annotation('textbox', [pos2(1)-0.09, pos2(2)+pos2(4)+0.1, 0.03, 0.03], ...
    'String','B','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

exportgraphics(gcf, fullfile(figureFolder,'F4ERPs_CrossCorrelation.png'),'Resolution', 600);
close
%% Figure 4: Mass-univariate results for TFRs, and cross correlation
% Load data from TFRs after averaging sessions from same subject
load(fullfile(dataFolder,filesep,'averagedTF.mat'))
% Load cross correlation results
load(fullfile(dataFolder,filesep,'cosineSimilarityPlotTF.mat'))
% Prepare
pairNames = {'blinkOffset_vs_fixation'   plot_blinkOff_fixation
    'blinkOffset_vs_blinkOnset' plot_blinkOff_blinkOn   
    'blinkOffset_vs_saccade'    plot_blinkOff_saccade
    'fixation_vs_blinkOnset'    plot_fixation_blinkOn
    'fixation_vs_saccade'       plot_fixation_saccade
    'blinkOnset_vs_saccade'     plot_blinkOn_saccade};
pairStr = pairNames(:,1);
pairVal = cell2mat(pairNames(:,2));
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
% Compute confidence interval (95%)
ciBlinkOff = compute_CI(finalBlinkOff,.05);
ciBlinkOn = compute_CI(finalBlinkOn,.05);
ciSacc = compute_CI(finalSaccade,.05);
ciFix = compute_CI(finalFixation,.05);
fontSize = 22;

figure(figure('Renderer', 'painters', 'Position', [10 10 1900 1200]))
% First, we plot the mass univariated results
t = tiledlayout(1,2);
t.TileSpacing = "compact";
t.Padding = "compact";
ax1 = nexttile;
hold (ax1,'on')
hold on
[~,~] = plotCI(squeeze(mean(finalSaccade(:,11,:),1)),squeeze(ciSacc(11,:,:)),time,cmapERP(35,:),fontSize,'-');
[~,~] = plotCI(squeeze(mean(finalBlinkOn(:,11,:),1)),squeeze(ciBlinkOn(11,:,:)),time,cmapERP(35,:),fontSize,'--');
[~,~] = plotCI(squeeze(mean(finalFixation(:,11,:),1)),squeeze(ciFix(11,:,:)),time,cmapERP(55,:),fontSize,'-');
[~,~] = plotCI(squeeze(mean(finalBlinkOff(:,11,:),1)),squeeze(ciBlinkOff(11,:,:)),time,cmapERP(55,:),fontSize,'--');
xlim([min(time) max(time)])
legend('Saccade Onset','','Blink Onset','','Fixation Onset','','Blink Offset')
v = vline(0,'-');
v.Color = [0 0 0 0.8];
v.LineWidth = 1.5;
h = hline(0,'-');
h.Color = [0 0 0 0.8];
h.LineWidth = 1.5;
xlabel('Time (s)')
ylabel('Power (dB)')
set(gca,'fontsize',fontSize)
% Now we plot the results from the cross correlation
ax2 = nexttile;
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
pos2 = ax2.Position;

annotation('textbox', [pos1(1)-0.03, pos1(2)+pos1(4), 0.03, 0.03], ...
    'String','A','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

annotation('textbox', [pos2(1)-0.09, pos2(2)+pos2(4)+0.1, 0.03, 0.03], ...
    'String','B','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

exportgraphics(gcf, fullfile(figureFolder,'F5TFRs_CrossCorrelation.png'),'Resolution', 600);
close
%% Figure 6: ERP locked to fixation onset and sorted by saccade duration
% Load the data
load(fullfile(rawDataFolder,filesep,'binnedData.mat'))
% Define parameters
subjects = [1,2,5,7,27,28];
titles = {'Subject 9 - Session 1','Subject 9 - Session 2',...
    'Subject 11 - Session 2','Subject 12 - Session 2',...
    'Subject 33 - Session 2','Subject 35 - Session 1'};
figure(figure('Renderer', 'painters', 'Position', [10 10 1000 800]))
t = tiledlayout(3,2,'TileSpacing', 'Compact', 'Padding', 'Compact');
for i = 1:length(subjects)
    nexttile(i)
    f = imagesc(times, 1:size(binnedDataALL{1,subjects(i)}, 2), ...
        binnedDataALL{1,subjects(i)}');
    colormap(magma(5000))
    xline(0, 'k-', 'LineWidth', 2);
    col = colorbar;
    col.Label.String = 'μV';
    maxVal = max(abs(binnedDataALL{1,subjects(i)}(:)));
    clim([-maxVal maxVal])
    xlabel('Time (ms)','FontSize',14);
    ylabel('Bins of Trials','FontSize',14);
    hold on;
    nBins = size(binnedDataALL{1,subjects(i)}, 2);
    scatter(saccadeOnsetsALL{1,subjects(i)}, 1:nBins, 14, 'w', 'filled');
    hold off;
    title(titles{i}, 'Interpreter', 'none');
    set(gca,'fontsize',14)
end

exportgraphics(gcf, fullfile(cd,'F6FixationTrialsSorted_Cz.png'), ...
    'Resolution', 600);
close