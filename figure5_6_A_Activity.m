%% Figure 5 and 6 A. Time-series (5) and Power (6)
%% Paths
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%% Figure 5A
% Load data from ERPs after averaging sessions from same subject
dataFolder = [projectFolder,filesep,'Data',filesep,'unfold'];
load(fullfile(dataFolder,filesep,'averagedERP.mat'))
% Define parameters
cmapERP = magma(70);
% Compute confidence interval (95%)
ciBlinkOff = compute_CI(finalBlinkOff,.05);
ciBlinkOn = compute_CI(finalBlinkOn,.05);
ciSacc = compute_CI(finalSaccade,.05);
ciFix = compute_CI(finalFixation,.05);
fontSize = 22;
figure
hold on
[~,~] = plotCI(squeeze(mean(finalSaccade(:,11,:),1)),squeeze(ciSacc(11,:,:)),time,cmapERP(35,:),fontSize,'-');
[~,~] = plotCI(squeeze(mean(finalBlinkOn(:,11,:),1)),squeeze(ciBlinkOn(11,:,:)),time,cmapERP(35,:),fontSize,'--');
[~,~] = plotCI(squeeze(mean(finalFixation(:,11,:),1)),squeeze(ciFix(11,:,:)),time,cmapERP(55,:),fontSize,'-');
[~,~] = plotCI(squeeze(mean(finalBlinkOff(:,11,:),1)),squeeze(ciBlinkOff(11,:,:)),time,cmapERP(55,:),fontSize,'--');
xlim([min(time) max(time)])
legend('Saccade Onset','','Blink Onset','','Fixation Onset','','Blink Offset','Box','off')
v = vline(0,'-');
v.Color = [0 0 0 0.8];
v.LineWidth = 1.5;
h = hline(0,'-');
h.Color = [0 0 0 0.8];
h.LineWidth = 1.5;
xlabel('Time [s]')
ylabel('Amplitude [µV]')
set(gca,'fontsize',fontSize)

exportgraphics(gcf, fullfile(figureFolder,'Figure5A.png'),'Resolution', 600);
close
%% Figure 6A
% Load data from TFRs after averaging sessions from same subject
load(fullfile(dataFolder,filesep,'averagedTF.mat'))
% Compute confidence interval (95%)
ciBlinkOff = compute_CI(finalBlinkOff,.05);
ciBlinkOn = compute_CI(finalBlinkOn,.05);
ciSacc = compute_CI(finalSaccade,.05);
ciFix = compute_CI(finalFixation,.05);
fontSize = 22;
figure
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
xlabel('Time [s]')
ylabel('Power [dB]')
set(gca,'fontsize',fontSize)
exportgraphics(gcf, fullfile(figureFolder,'Figure6A.png'),'Resolution', 600);
close