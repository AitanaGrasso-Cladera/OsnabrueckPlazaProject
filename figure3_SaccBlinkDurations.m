%% Figure 3. Averaged distribution of blink and saccade durations
%%
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%% For every session, extract the duration of each valid saccade and blink
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
%% Plot into histograms
colors = viridis(50);
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
    data = saccDur{sub};
    data(data == 0) = [];
    % histogram per subject (normalized)
    counts = histcounts(data, binEdges, 'Normalization', 'probability');
    counts = smoothdata(counts, 'gaussian', 8);
    allCounts(sub,:) = counts;
end
% Mean histogram
meanCounts = mean(allCounts, 1);
plot(binCenters, meanCounts, ...
    'k', ...
    'LineWidth', 2);
xlabel('Duration [ms]');
ylabel('Frequency of occurence');
grid on;
hold off;
a = gca;
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
ylim([0 0.13])
set(gca,'fontsize',14)
xtk = get(gca, 'XTick'); 
set(gca, 'XTick', (0:100:500))
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
    % bar(binCenters, counts, ...
    %     'FaceAlpha', 0.12, ...
    %     'EdgeColor', 'none', ...
    %     'FaceColor', colors(30,:));
end

% Mean histogram
meanCounts = mean(allCounts, 1);
plot(binCenters, meanCounts, ...
    'k', ...
    'LineWidth', 2);
xlabel('Duration [ms]');
ylabel('Frequency of occurence');
grid on;
hold off;
a = gca;
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
ylim([0 0.13])
set(gca,'fontsize',14)
xtk = get(gca, 'XTick'); 
set(gca, 'XTick', (0:100:500))

pos1 = ax1.Position;
pos2 = ax2.Position;
fontSize = 16;
annotation('textbox', [pos1(1)-0.07, pos1(2)+pos1(4)+0.02, 0.03, 0.03], ...
    'String','A','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')
annotation('textbox', [pos2(1)-0.07, pos2(2)+pos2(4)+0.02, 0.03, 0.03], ...
    'String','B','LineStyle','none', ...
    'FontSize',fontSize+2,'FontWeight','bold')

exportgraphics(gcf, fullfile(figureFolder,'Figure3.png'),'Resolution', 600);
close