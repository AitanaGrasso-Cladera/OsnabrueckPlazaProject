%% Compute correlation betweeen theoretical matrix and each result from the similarity 
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
    load([dataFolder,filesep,'finalDatasetCosineERP.mat'])
elseif POWER == 1
    load([dataFolder,filesep,'finalDatasetCosineTF.mat'])
end
%% Generate a matrix per subject and channel
for sub = 1:size(final_fixation_saccade,1)
    for i = 1:size(final_fixation_saccade,2)
        matrices{sub,i} = [1,final_blinkOff_fixation(sub,i),final_blinkOff_blinkOn(sub,i),final_blinkOff_saccade(sub,i);...
            final_blinkOff_fixation(sub,i),1,final_fixation_blinkOn(sub,i),final_fixation_saccade(sub,i);...
            final_blinkOff_blinkOn(sub,i),final_fixation_blinkOn(sub,i),1,final_blinkOn_saccade(sub,i);...
            final_blinkOff_saccade(sub,i),final_fixation_saccade(sub,i),final_blinkOn_saccade(sub,i),1];
    end
end
T = [1,1,-1,-1;1,1,-1,-1;-1,-1,1,1;-1,-1,1,1];
template = T;  
ut = triu(true(size(template)),1);  % mask for upper triangle
template_vec = template(ut);

similarity = zeros(16,11);
for sub = 1:16
    for i = 1:11
        M = matrices{sub,i};  % your 11 matrices
        similarity(sub,i) = corr(M(ut), template_vec, 'type', 'Pearson');
    end
end
%% Plot with topographical distribution
% Electrode positions
x = [-0.3 0 0.3 -0.57 -0.5  -0.5     0.5   0.5  0.57  -0.57 0.57];
y = [  0  0  0  -0.2  -0.1  0.05    0.05  -0.1  -0.2    0.2  0.2];

vals = mean(similarity,1);

figure
hold on

% Draw head outline
th = linspace(0,2*pi,300);
plot(0.6*cos(th),0.6*sin(th), ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',2);

% Draw nose
plot([-0.05 0 0.05], ...
     [0.60 0.66 0.60], ...
     'Color',[0.4 0.4 0.4], ...
     'LineWidth',2);

% Draw ears
t = linspace(-0.8*pi/2,0.8*pi/2,100);

plot(-0.60-0.03*cos(t),0.11*sin(t), ...
     'Color',[0.4 0.4 0.4],'LineWidth',2)

plot( 0.60+0.03*cos(t),0.11*sin(t), ...
     'Color',[0.4 0.4 0.4],'LineWidth',2)

% Electrodes
scatter(x,y,400,vals,'filled');

colormap(magma)
clim([0 1])

c = colorbar;
c.Label.String = 'Pearsons correlation (r)';
c.Ticks = [0 1];
c.FontSize = 16;

labels = {'C3','Cz','C4','L4','L3','L2','R2','R3','R4','L1','R1'};

for k = 1:length(labels)
    text(x(k),y(k)+0.06,labels{k}, ...
         'HorizontalAlignment','center', ...
         'FontWeight','bold');
end

axis equal
axis off

xlim([-0.8 0.8])
ylim([-0.75 0.75])

hold off
%% Save
if ERP == 1
    save([dataFolder,filesep,'correlationERP.mat'],'similarity')
    exportgraphics(gcf, fullfile(cd,'Figure5E.png'),'Resolution', 600);
elseif POWER == 1
    save([dataFolder,filesep,'correlationTF.mat'],'similarity')
    exportgraphics(gcf, fullfile(cd,'Figure6E.png'),'Resolution', 600);
end
close
