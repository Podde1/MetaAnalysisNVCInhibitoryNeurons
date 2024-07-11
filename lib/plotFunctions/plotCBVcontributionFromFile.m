function [] = plotCBVcontributionFromFile()

load('Parameters/CellSpecificContributions/CBVTable.mat', 'CBVTable')

%% Load and sort data
Experiment   = cell( size(CBVTable{1,2},1)*size(CBVTable(1,2:end),2), 1, size(CBVTable,1) );
NeuronType   = cell( size(CBVTable{1,2},1)*size(CBVTable(1,2:end),2), 1, size(CBVTable,1) );
Contribution = zeros(size(CBVTable{1,2},1)*size(CBVTable(1,2:end),2), 1, size(CBVTable,1) );
NeuronTypes = {'NO', 'NPY', 'Pyr', 'Ast'};
ShortExperimentNames = {'Vazques Sens1 4s', 'Vazques Exc 5Hz 4s', 'Vazques Sens2 4s', 'Vazques Inh 5Hz 4s', 'Moon Exc 1Hz 20s', 'Moon Exc 20Hz 20s',  'Moon Sens 20s', 'Moon Inh 1Hz 20s', 'Moon Inh 20Hz 20s','Moon Inh 1Hz 5s', 'Moon Inh 20Hz 5s'};
SimpleTitles = {'Experiment 1', 'Experiment 2', 'Experiment 3', 'Experiment 4', 'Experiment 5', 'Experiment 6', 'Experiment 7', 'Experiment 8', 'Experiment 9', 'Experiment 10', 'Experiment 11'};

%sort table by stimuli paradigm
rearrangeOrder = [7 1 3 2 5 6 4 8 10 9 11];
CBVTable = CBVTable(rearrangeOrder,:);
ShortExperimentNames = ShortExperimentNames(:,rearrangeOrder);

for i = 1:size(CBVTable,1) % rows 
    count = 1;
    for j = 1:size(CBVTable(1,2:end),2) % columns
        for k = 1:size(CBVTable{1,2},1) % elements in array
            Experiment{count,1,i}   = CBVTable{i};
            NeuronType{count,1,i}   = NeuronTypes{j};
            Contribution(count,1,i) = CBVTable{i,j+1}{k};
            count = count + 1;
        end
    end
end

ColorAst = [0 168 0]./255;
ColorNO  = [237 177 32]./255; 
ColorNPY = [217 83 25]./255;
ColorPyr = [87 155 255]./255;
Neuroncolors = [ColorNO; ColorNPY; ColorPyr; ColorAst];

ExpColorExc  = [174 208 255 36*2.55]./255;
ExpColorSens = [229 229 229 100*2.55]./255;
ExpColorInhLow  = [239 153 119 36*2.55]./255;
ExpColorInhHigh  = [241 199 97 36*2.55]./255;
ExperimentColors = [ExpColorSens; ExpColorExc; ExpColorSens; ExpColorInhLow; ExpColorExc; ExpColorExc; ExpColorSens; ExpColorInhLow; ExpColorInhHigh; ExpColorInhLow; ExpColorInhHigh];

%rearrange the colour order 
ExperimentColors = ExperimentColors(rearrangeOrder,:);

numElements = size(ExperimentColors,1);

%% Boxplot
nrows = 2;
ncols = 6; 

figure()
dpi = 96;
a4WidthPixels = round(210 / 25.4 * dpi);
a4HeightPixels = round(0.4*297 / 25.4 * dpi);
set(gcf, 'Position', [0 0 a4WidthPixels a4HeightPixels]);

warning('off')
tiledlayout(nrows, ncols, 'Padding', 'none', 'TileSpacing', 'compact')
for i = 1:numElements
    
    nexttile
    bp = boxplot(Contribution(:,1,i), NeuronType(:,1,i), 'Symbol', 'o', 'GroupOrder', NeuronTypes, 'Colors', Neuroncolors, 'Whisker', 10);
   
    ax = gca;
    ax.Color = ExperimentColors(i,1:4);

    % Get handles to the box plot elements
    h = findobj(bp, '-property', 'Tag');
    
    % Update objects
    set(ax,{'YLim','Color','XColor','YColor','LineWidth','FontSize', 'Box'},{[0 100], ExperimentColors(i,1:4), [0 0 0], [0 0 0], 2, 16, 'off'})
    
    % set properties of Whisker caps.
    set(h(or(strcmp(get(h,'Tag'),'Upper Adjacent Value'),strcmp(get(h,'Tag'),'Lower Adjacent Value'))),{'Color'},{[0 0 0]})
    hit = 0;
    for j = 1:length(h)
        if or(strcmp(get(h(j), 'Tag'), 'Upper Whisker'), strcmp(get(h(j), 'Tag'), 'Lower Whisker')) 
            % set line style
            set(h(j), {'LineStyle','Color'}, {'-', [0 0 0]});

        elseif strcmp(get(h(j), 'Tag'), 'Box')
            % create a rectangle to fill in the box
            hit = hit + 1;
            recWidth  = h(j).XData(3)- h(j).XData(1);
            recHieght = h(j).YData(3)- h(j).YData(1);

            rectangle(ax, 'Position', [h(j).XData(1) h(j).YData(1) recWidth recHieght],...
              'FaceColor', Neuroncolors(hit, 1:3), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        end

        % Set line widths 
        set(h(j), 'LineWidth', 1.5);
    end

    title(SimpleTitles{i}, 'Interpreter', 'none', 'FontSize', 12)
    ylim([0 100])
    ylabel('Contribution (%)', 'FontSize', 4)
    yticks([0 20 40 60 80 100])
    yticklabels([0 20 40 60 80 100])
    set(ax.XAxis, 'FontSize', 10)
    set(ax.YAxis, 'FontSize', 10)
    xtickangle(45)
end

% Crete legends
nexttile
xticks([])
yticks([])

for i = 1:numElements
    line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color', ExperimentColors(i,:), 'DisplayName', ShortExperimentNames{i}) %line to be able to rectangle add to legend
end

rectHandles = findobj(gca, '-property', 'Tag', 'Type', 'Line');
idxExc     = 0;
idxInhLow  = 0; 
idxInhHigh = 0; 
idxSens    = 0;

for idxRectHandle = 1:length(rectHandles)
    if strcmp(rectHandles(idxRectHandle).DisplayName, 'Moon Exc 1Hz 20s')
        idxExc = idxRectHandle;
    elseif strcmp(rectHandles(idxRectHandle).DisplayName, 'Moon Sens 20s')
        idxSens = idxRectHandle;
    elseif strcmp(rectHandles(idxRectHandle).DisplayName, 'Moon Inh 1Hz 20s')
        idxInhLow = idxRectHandle;
    elseif strcmp(rectHandles(idxRectHandle).DisplayName, 'Moon Inh 20Hz 20s')
        idxInhHigh = idxRectHandle;
    end
end

[~, l] = legend(rectHandles([idxExc, idxInhLow, idxInhHigh, idxSens]),{'Excitatory', 'Inhibitory Low', 'Inhibitory Low', 'Sensory'}, 'FontSize', 12);
objhl = findobj(l, 'type', 'line'); 
set(objhl, 'Linewidth', 8); 

exportgraphics(gcf,'Figures/BoxplotCBV.png', 'Resolution', 600)
exportgraphics(gcf,'Figures/BoxplotCBV.pdf', 'ContentType', 'vector')


%% Spiderplot
figure()
a4HeightPixels = round(0.6*297 / 25.4 * dpi);
set(gcf, 'Position', [0 0 a4WidthPixels a4HeightPixels]);

spiderMatrix = zeros(size(CBVTable,1), size(CBVTable(1,2:end),2));

for i=1:size(spiderMatrix,1)
    for j=1:size(spiderMatrix,2)
        spiderMatrix(i,j) = CBVTable{i,j+1}{1};
    end
end

spider_plot(spiderMatrix', 'axeslimits', repmat([0 100], 11, 1)', 'Color', Neuroncolors, 'axeslabels', SimpleTitles, 'axeslabelsedge', 'none', 'filloption', 'on', ...
    'axesinterval', 1, 'axesfontsize', 12, 'labelfontsize', 16, 'axesoffset', 0, 'axesprecision', 0, 'axestickformat', '%.0f %%', ... 
    'BackgroundColor', [1 1 1], 'axesfontcolor', [0 0 0], 'labelcolor', [0 0 0], 'axeslabelsrotate', 'off');

for i = 1:4
    line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color', Neuroncolors(i,:), 'DisplayName', NeuronTypes{i}) %line to be able to rectangle add to legend
end

LineHandles = findobj(gca, '-property', 'Tag', 'Type', 'Line');
[~, l] = legend(LineHandles([4 3 2 1]), 'TextColor', [0 0 0], 'Edgecolor', 'none', 'FontSize', 12);
objhl = findobj(l, 'type', 'line'); 
set(objhl, 'Linewidth', 8);

warning('on')

exportgraphics(gcf,'Figures/SpiderCBV.png', 'Resolution', 600)
exportgraphics(gcf,'Figures/SpiderCBV.pdf', 'ContentType', 'vector')

end