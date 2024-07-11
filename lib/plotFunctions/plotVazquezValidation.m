function [] = plotVazquezValidation()

load('Parameters/Vazquez_2ms_10ms_Uncertainty.mat', 'simulationUncertainty')

modelName = 'Sten2023_v21';
ExperimentIdx = [2 3 6 7];
VarNames = {'HbT', 'HbR', 'HbO', 'CBF', 'BOLD_OIS'};
[~, DataStructure, ~, variableNames, ~, ~, ~, ~] = modelSetup(modelName,ExperimentIdx,VarNames); 

%% Plotting
simulationFieldNames = fieldnames(simulationUncertainty); % get the field names of the simulation structure. These will correspond to the simulated experiment field names
nVariables = cellfun('size',variableNames,1); % returns a vectro with the number of variables for each experiment 

%% Set up figures
Vazquez_HbX_figure = figure('Units','centimeters','Position',[0 0 21.0 29.7]);
Vazquez_HbX_ax = tiledlayout(4, 3, "TileSpacing","tight", "Padding","tight"); % preallocate varaible for axis object for Vazquez plot
Vazquez_HbX_ax_Titles = {'HbO', 'HbR', 'HbT', 'CBF', 'BOLD-OIS', '', 'HbO', 'HbR', 'HbT', 'CBF', 'BOLD-OIS', ''};

globalFontSize = 12;

for i = [1:5 7:11]
    Vazquez_HbX_ax(i) = nexttile(i);

    set(Vazquez_HbX_ax(i),{'XLim','Fontsize'},{[-10 60],globalFontSize})
    set(Vazquez_HbX_ax(i).Title,{'String'},{Vazquez_HbX_ax_Titles{i}})
    set(Vazquez_HbX_ax(i).YLabel,{'String'},{['\Delta ',Vazquez_HbX_ax_Titles{i} ,'(%)']})
    set(Vazquez_HbX_ax(i).XLabel,{'String'},{'Time (s)'})
end

PlotStructure.Vazquez2018_HbX = {'Vazques2018_Excitatory_FL'  ...
                                 'Vazques2018_Excitatory_2ms' ...
                                 'Vazques2018_Excitatory_10ms'...
                                 'Vazques2018_Excitatory_30ms'...
                                 'Vazques2018_Inhibitory_FL'  ...
                                 'Vazques2018_Inhibitory_2ms' ...
                                 'Vazques2018_Inhibitory_10ms'...
                                 'Vazques2018_Inhibitory_30ms' };

Colour = {[0.3     0.3     0.3         % Vazques2018_Excitatory_FL         
           0.5020  0.6980  1           % Vazques2018_Excitatory_2ms
           0       0.4470  0.7410      % Vazques2018_Excitatory_10ms
           0       0.2627  0.6510      % Vazques2018_Excitatory_30ms
           0.5     0.5     0.5         % Vazques2018_Inhibitory_FL
           1       0.8     0.5020      % Vazques2018_Inhibitory_2ms
           0.9290  0.6940  0.1250      % Vazques2018_Inhibitory_10ms
           0.8500  0.3250  0.0980]};   % Vazques2018_Inhibitory_30ms
       

plotVariableNames = {strrep(Vazquez_HbX_ax_Titles([1 2 3 4 5 7 8 9 10 11]),'-','_')};
plotIdx           = {[repmat([intersect(ExperimentIdx,[2 6])],5,1); repmat([intersect(ExperimentIdx,[3 7])],5,1)]}; 
                                
plotStructureFieldNames = fieldnames(PlotStructure);
figureHandles = [Vazquez_HbX_figure];

% plot the data
errorBarHandle = mat2cell(gobjects(sum(cellfun(@numel,plotIdx,'uni',1)),1),cellfun(@numel,plotIdx,'uni',1),[1]);
errorBarHandle = cellfun(@reshape,errorBarHandle,cellfun(@(x)fliplr(size(x)),plotIdx,'uni',0).','uni',0);
for i = 1:numel(figureHandles)% loop over the three main figures
    currentAxes = flipud(findobj(get(figureHandles(i),'Children'), '-depth', 1, 'type', 'axes')); % get all axis objects from figureHandel i 
    for j = 1:numel(currentAxes) % loop over the mumber of axes in figure i
        line(currentAxes(j),currentAxes(j).XLim,[0 0],'LineWidth',0.5,'Color',[0.8 0.8 0.8])% add a thin line along y=0 to indicate baseline
        count_k = 1;
        for k = plotIdx{i}(j,:) %numel(PlotStructure.(plotFieldNames{i})) % loop over the by plotIdx indicated experiments of PlotStructure field i.
            currentDataStructure = DataStructure.(PlotStructure.(plotStructureFieldNames{i}){k}).(plotVariableNames{i}{j}); % get the data strucute that corresponts to plotStructure: field i, entry k, varaible j. 
            hold(currentAxes(j),'on')
            errorBarHandle{i}(count_k,j) = errorbar(currentAxes(j),currentDataStructure.Time, currentDataStructure.Mean, currentDataStructure.SEM, ' *', 'Color',Colour{i}(k,:),'LineWidth',1, 'CapSize', 2, 'MarkerSize', 4);
            count_k = count_k + 1;
        end
    end
end

%% plot the simulations
plotsimulationUncertainty = 1;

for i = 1:numel(simulationFieldNames) % loop over the simualtion structure fields i.e. the experiments we have simualtions for
    for j = 1:numel(plotStructureFieldNames) % loop over the differnt plot structure feilds 
        % Find the index of simulation field i in plof structure field j. 
        simulationFieldInPlotFieldIdx = find(strcmp(simulationFieldNames{i},PlotStructure.(plotStructureFieldNames{j})),1);
        if ~isempty(simulationFieldInPlotFieldIdx) % if simulation field is in plot field j 
            % Get the axes objects of plot j. 
            currentAxes = flipud(findobj(get(figureHandles(j),'Children'), '-depth', 1, 'type', 'axes'));
            for k = 1:nVariables(i) % loop over the differnt varaibles in simulation field i.
                % get the logical index of variable k in plotIdex cell j. Check if varaible name k is in plotVaraibleName cell j, and if simulation field i should be plotted in plot j (relevant for Moon Hbx and Desjardins). 
                variableIdx = and(strcmp(variableNames{i}(k),plotVariableNames{j}),any(plotIdx{j}==simulationFieldInPlotFieldIdx,2).');
                if any(variableIdx) % if the logical index has any ture values plot the simualtion i, variable k in plot j 
                    if plotsimulationUncertainty
                        fill(currentAxes(variableIdx),[simulationUncertainty.(simulationFieldNames{i}).t,fliplr(simulationUncertainty.(simulationFieldNames{i}).t)],[simulationUncertainty.(simulationFieldNames{i}).(variableNames{i}{k}).simulation_lb,fliplr(simulationUncertainty.(simulationFieldNames{i}).(variableNames{i}{k}).simulation_ub)],'',...
                            'FaceColor',Colour{j}(simulationFieldInPlotFieldIdx,:),...
                            'FaceAlpha',0.2, ...
                            'EdgeColor',Colour{j}(simulationFieldInPlotFieldIdx,:))
                    end
                end
            end
        end
    end
end

%% set legends
legendStrings ={{'Excitatory FL'  ;
                'Excitatory 2ms' ;
                'Excitatory 10ms';
                'Excitatory 30ms';
                'Inhibitory FL'  ;
                'Inhibitory 2ms' ;
                'Inhibitory 10ms';
                'Inhibitory 30ms'}};

for i = 1:numel(figureHandles)
    currentAxes = findobj(get(figureHandles(i),'Children'), '-depth', 1, 'type', 'axes');
    count_k = 0;
    for k = ExperimentIdx
        count_k = count_k + 1;
        legLineHandle{i}(count_k) = line(currentAxes(1),[-11 -10.5], [0 0],'LineWidth',5,'Color',[Colour{i}(k,:)]);
    end
    legend(legLineHandle{i},legendStrings{i}(ExperimentIdx),'AutoUpdate','off','FontSize',globalFontSize)
end

%% plot the line for the stimulation periond after everything else to align with the x-axis
for i = 1:numel(figureHandles) % loop over the three main figures
    currentAxes = flipud(findobj(get(figureHandles(i),'Children'), '-depth', 1, 'type', 'axes')); % get all axis objects from figureHandel i 
    for j = 1:numel(currentAxes) % loop over the mumber of axes in figure i
        k = plotIdx{i}(j,end);
        line(currentAxes(j),[0 DataStructure.(PlotStructure.(plotStructureFieldNames{i}){k}).Stimulation.Length],repmat(currentAxes(j).YLim(1),1,2),'LineWidth',5,'Color',[0 0 0])% add a thin line at thebottomof the plo to illustrate stimulation time
    end
end

%% add A-J subfigure notations to the figures 
addSubFigureNotations = 1;
figLabelFontsize = 16;

if addSubFigureNotations
    setSubFigureNotations(Vazquez_HbX_ax([1:5 7:11]), figLabelFontsize);
end

end