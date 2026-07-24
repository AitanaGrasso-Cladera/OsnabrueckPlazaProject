%% Figure 5 and 6 D. Theoretical matrix
%% Paths
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%%
fontSize = 16;
imagesc(T)
hold on
axis equal tight
set(gca,'XTick',1:n,'XTickLabel',displayLabels,'XTickLabelRotation',45,...
        'YTick',1:n,'YTickLabel',displayLabels,'YDir','reverse');
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

textStrings = num2str(T(:),'%0.2f');
textStrings = strtrim(cellstr(textStrings));

[x,y] = meshgrid(1:n);

hText = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center',...
    'FontSize',fontSize,'FontWeight','bold');

midpoint = 0;

for k = 1:numel(T)
    if T(k) > midpoint
        set(hText(k),'Color','k')
    else
        set(hText(k),'Color','w')
    end
end

hold on
xline(2.5,'k','LineWidth',2)
yline(2.5,'k','LineWidth',2)
set(gca,'fontsize',fontSize)

exportgraphics(gcf, fullfile(figureFolder,'Figure5_6_theoreticalMatrix.png'),'Resolution', 600);
close