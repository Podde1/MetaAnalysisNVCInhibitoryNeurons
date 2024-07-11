function [figureHandles] = plotCurrenSolution(theta, modelName,ExperimentIdx,variables2Fit,plotIndividualExperiments,plotAllStates)
% close all
if nargin < 5; plotIndividualExperiments = 0; end
if nargin < 6; plotAllStates=0;end

addpath(genpath('./Models'))

%% run model setup 
[circuitConstants, DataStructure, model, variableNames, ~, modelObservablesIndex, parameterMappingTable, options] = modelSetup(modelName,ExperimentIdx,variables2Fit);

%% simulate model and calculate fval
thetaConstants = []; OptParamBool=[]; %set to empty since we will always wannna use the entire parameter vector for the ploting 

[fval,simulation] = objectiveFunction(theta,modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure, ExperimentIdx, variableNames, modelObservablesIndex, parameterMappingTable, model, options,0,[]);

fprintf('Objective function value for current solution:\t%0.4f\n',fval)

experimentFieldNames = fieldnames(DataStructure);

%% plot all fitted variables as subplots with one figure for each simulation/experiment
simulationFieldNames = fieldnames(simulation); % get the field names of the simulation structure. These will correspond to the simulated experiment field names
nVariables = cellfun('size',variableNames,1); % returns a vector with the number of variables for each experiment 
ax = mat2cell(gobjects(sum(nVariables),1),nVariables,1); % praallocates the cell array of axes-handles as generic graphical object
ylableStr = strcat('\Delta ',strrep(variableNames{1},'_','-'),' (%)');
count = 1;
if plotIndividualExperiments
    for i = 1:numel(simulationFieldNames) % loop over all simulations
        f1=figure(); % create a new figure for the i:th simulation

        % calculate the squarest numeber of subplots needed
        numberOfPlots = numel(variableNames{i});
        nRows = floor(sqrt(numberOfPlots));
        nCols =  ceil(sqrt(numberOfPlots));
        if nRows*nCols < numberOfPlots
            nCols = nCols + 1;
        end

        for j=1:numel(variableNames{i}) % loop over all fitted variables in the i:th simulation
            ax{i}(j) = subplot(nRows,nCols,j,'Parent',f1);
            hold(ax{i}(j),'on')
            TempDataStruct = DataStructure.(experimentFieldNames{ExperimentIdx(i)}).(variableNames{i}{j}); % create a temporary data structure to avoid repetitive calling of the same structure.
            erbarHandle = errorbar(TempDataStruct.Time,TempDataStruct.Mean,TempDataStruct.SEM,' r*','LineWidth',2, 'CapSize', 2, 'MarkerSize', 4);  % plot data as red error bars
            simplotHandle = plot(simulation.(simulationFieldNames{i}).t,simulation.(simulationFieldNames{i}).y(:,modelObservablesIndex{count}(j)),'b-','LineWidth',2); % plot simulation as blue line

            line([0, 4],repmat(ax{i}(j).YLim(1),1,2),'LineWidth',5,'Color',[0 0 0]) % plot a black line at the bottom of the axes to representing the stimulation time.
            title(strrep(strjoin(variableNames{i}(j),' '),'_',' ')) % set the title of axes j, to the corresponding variable name.
            ylabel(ylableStr{j});
            xlabel('Time (s)');
            hold(ax{i}(j),'off')
        end
        sgtitle(f1,strrep(experimentFieldNames{ExperimentIdx(i)},'_',' ')) % set the sub-group titel to the experiment name

        count = count + 1;
    end
end

if plotAllStates
    for i = 1:numel(simulationFieldNames) % loop over all simulations
        if and(nargin >=5, plotAllStates)
            % plot all states for each simulation in a seperate figure
            fAllStates = plotStates(simulation.(simulationFieldNames{i}),strrep([model.states, model.observables],'_','-'));
            sgtitle(fAllStates,strrep(simulationFieldNames{i},'_',' '))
        end
    end
end
%% Set up figures
plotFullScreen = 1;

if plotFullScreen; windowStateStr = 'maximize';else; windowStateStr = 'normal';end

BOLD_fMRI_figure = figure('Units','normalized','WindowState',windowStateStr);
BOLD_fMRI_ax = axes(BOLD_fMRI_figure);
set(BOLD_fMRI_ax,{'XLim','YLim','Fontsize'},{[-10 72],[-3.5 6],24})
set(BOLD_fMRI_ax.Title,{'String'},{'Moon 2021 BOLD-fMRI'})
set(BOLD_fMRI_ax.YLabel,{'String'},{'\Delta BOLD-fMRI (%)'})
set(BOLD_fMRI_ax.XLabel,{'String'},{'Time (s)'})

Moon_HbX_figure = figure('Units','normalized','WindowState',windowStateStr);
Vazquez_HbX_figure = figure('Units','normalized','WindowState',windowStateStr);
Lee_HbX_figure = figure('Units','normalized','WindowState',windowStateStr);

Moon_HbX_ax = gobjects(6,1); % preallocate varaible for axis object for Monn HbX plot
Vazquez_HbX_ax = gobjects(6*2,1); % preallocate varaible for axis object for Vazquez plot
Lee_HbX_ax = gobjects(9,1); % preallocate varaible for axis object for Lee plot

Vazquez_HbX_insertAx = gobjects(5,1);
Lee_HbX_insertAx = gobjects(3,1);


Moon_HbX_ax_Titles = {'HbO', 'HbR', 'HbT', 'HbO', 'HbR', 'HbT'};
Vazquez_HbX_ax_Titles = {'HbO', 'HbR', 'HbT', 'CBF', 'BOLD-OIS', ''};
Lee_HbX_ax_Titles = {'HbO', 'HbR', 'HbT', 'HbO', 'HbR', 'HbT'};

globalFontSize = 16;
insertAxisScale = 0.5;
insertFontSize = 12;

for i = 1:6
    Moon_HbX_ax(i) = subplot(2,3,i,'Parent',Moon_HbX_figure);
    Vazquez_HbX_ax(i) = subplot(4,3,i+(ceil(i/3)-1)*3,'Parent',Vazquez_HbX_figure);
    if i<6; Vazquez_HbX_insertAx(i) = subplot(4,3,(i+3)+(ceil(i/3)-1)*3,'Parent',Vazquez_HbX_figure); end
    Lee_HbX_ax(i) = subplot(3,3,i,'Parent',Lee_HbX_figure);
     
    set(Moon_HbX_ax(i),{'XLim','Fontsize'},{[-7 62 - round(heaviside(i-4))*35],globalFontSize})
    set(Moon_HbX_ax(i).Title,{'String'},{Moon_HbX_ax_Titles{i}})
    set(Moon_HbX_ax(i).YLabel,{'String'},{['\Delta ',Moon_HbX_ax_Titles{i} ,'Conc (\muM)']})
    set(Moon_HbX_ax(i).XLabel,{'String'},{'Time (s)'})
    
    if i<6
        set(Vazquez_HbX_ax(i),{'XLim','Fontsize'},{[-10 60],globalFontSize})
        set(Vazquez_HbX_ax(i).Title,{'String'},{Vazquez_HbX_ax_Titles{i}})
        set(Vazquez_HbX_ax(i).YLabel,{'String'},{['\Delta ',Vazquez_HbX_ax_Titles{i} ,'(%)']})
        set(Vazquez_HbX_ax(i).XLabel,{'String'},{'Time (s)'})
        set(Vazquez_HbX_insertAx(i),{'XLim','XTickLabel','YTickLabel','FontSize'},{[-10 60],{},{},insertFontSize})
    end
    
    Vazquez_HbX_ax(i).Position(2) = Vazquez_HbX_ax(i).Position(2) - Vazquez_HbX_ax(i).Position(4);  % move the y position of the axis down by the hight of the axis  
    Vazquez_HbX_ax(i).Position(4) = Vazquez_HbX_ax(i).Position(4)*2; % double the height of the axis
    if i<6
        Vazquez_HbX_insertAx(i).Position([3,4]) = Vazquez_HbX_insertAx(i).Position([3,4]).*insertAxisScale; %shrink the insert axis 
        Vazquez_HbX_insertAx(i).Position([1,2]) = Vazquez_HbX_ax(i).Position([1,2]) + Vazquez_HbX_ax(i).Position([3,4]) - Vazquez_HbX_insertAx(i).Position([3,4]);%position the insert axis in the toop left corner of the big axis 
    end

    set(Lee_HbX_ax(i),{'XLim','Fontsize'},{[[-12 63] - round(heaviside(i-4))*[-9 49]],globalFontSize})
    set(Lee_HbX_ax(i).Title,{'String'},{Lee_HbX_ax_Titles{i}})
    set(Lee_HbX_ax(i).YLabel,{'String'},{['\Delta ',Lee_HbX_ax_Titles{i} ,'Conc (\muM)']})
    set(Lee_HbX_ax(i).XLabel,{'String'},{'Time (s)'})
end

for i = 1:3
    Lee_HbX_insertAx(i) = subplot(3,3,i+6,'Parent',Lee_HbX_figure);
    set(Lee_HbX_insertAx(i),{'XLim','XTickLabel','YTickLabel','FontSize'},{[-12 63],{},{},insertFontSize})
    Lee_HbX_insertAx(i).Position([3,4]) = Lee_HbX_insertAx(i).Position([3,4]).*insertAxisScale; %shrink the insert axis 
    Lee_HbX_insertAx(i).Position([1,2]) = Lee_HbX_ax(i).Position([1,2]) + Lee_HbX_ax(i).Position([3,4]) - Lee_HbX_insertAx(i).Position([3,4]); %position the insert axis in the toop left corner of the big axis 
end

for i = 1:6 % has to be done after all axes are created for the lee axes not to overwrite each other 
    Lee_HbX_ax(i).Position= Vazquez_HbX_ax(i).Position;
end

%%
PlotStructure.Moon2021_BOLD_fMRI = {'Moon2021_Excitatory_1Hz_20s' ...
                                    'Moon2021_Excitatory_20Hz_20s'...
                                    'Moon2021_Forepaw_4Hz_20s'    ...
                                    'Moon2021_Inhibitory_1Hz_20s' ...
                                    'Moon2021_Inhibitory_20Hz_20s'};

PlotStructure.Vazquez2018_HbX = {'Vazques2018_Excitatory_FL'  ...
                                 'Vazques2018_Excitatory_2ms' ...
                                 'Vazques2018_Excitatory_10ms'...
                                 'Vazques2018_Excitatory_30ms'...
                                 'Vazques2018_Inhibitory_FL'  ...
                                 'Vazques2018_Inhibitory_2ms' ...
                                 'Vazques2018_Inhibitory_10ms'...
                                 'Vazques2018_Inhibitory_30ms' };

PlotStructure.Moon2021_HbX = {'Moon2021_Inhibitory_1Hz_20s' ...
                              'Moon2021_Inhibitory_20Hz_20s'...
                              'Moon2021_Inhibitory_1Hz_5s'  ...
                              'Moon2021_Inhibitory_20Hz_5s' };

PlotStructure.Lee2020_HbX = {'Lee2020_Whiskers_5Hz_16s'...
                             'Lee2020_Whiskers_5Hz_2s'...
                             'Lee2020_SST_20Hz_16s'...
                             'Lee2020_SST_20Hz_2s'...
                             'Lee2020_nNOS_20Hz_16s'...
                             'Lee2020_nNOS_20Hz_2s'};

Colour = {[0       0.4470  0.7410      % Moon2021_Excitatory_1Hz_20s                          
           0       0.2627  0.6510      % Moon2021_Excitatory_20Hz_20s
           0.5     0.5     0.5         % Moon2021_Forepaw_4Hz_20s
           0.8500  0.3250  0.0980      % Moon2021_Inhibitory_1Hz_20s
           0.9290  0.6940  0.1250]...  % Moon2021_Inhibitory_20Hz_20s

          [0.3     0.3     0.3         % Vazques2018_Excitatory_FL         
           0.5020  0.6980  1           % Vazques2018_Excitatory_2ms
           0       0.4470  0.7410      % Vazques2018_Excitatory_10ms
           0       0.2627  0.6510      % Vazques2018_Excitatory_30ms
           0.5     0.5     0.5         % Vazques2018_Inhibitory_FL
           1       0.8     0.5020      % Vazques2018_Inhibitory_2ms
           0.9290  0.6940  0.1250      % Vazques2018_Inhibitory_10ms
           0.8500  0.3250  0.0980]...  % Vazques2018_Inhibitory_30ms

          [0.8500  0.3250  0.0980      % Moon2021_Inhibitory_1Hz_20s
           0.9290, 0.6940, 0.1250      % Moon2021_Inhibitory_20Hz_20s
           0.8500  0.3250  0.0980      % Moon2021_Inhibitory_1Hz_5s
           0.9290, 0.6940, 0.1250 ]... % Moon2021_Inhibitory_20Hz_5s
           
          [0.5     0.5     0.5         % Lee 2020_Whiskers 2s 
           0.5     0.5     0.5         % Lee 2020_Whiskers 16s
           0.9     0.2     0.08        % Lee 2020_SST 2s
           0.9     0.2     0.08        % Lee 2020_SST 16s
           0.9290, 0.6940, 0.1250      % Lee 2020 nNOS 2s
           0.9290, 0.6940, 0.1250 ]};  % Lee 2020 nNOS 16s1      


plotVariableNames = {{'BOLD_fMRI'},...
                     strrep(Vazquez_HbX_ax_Titles([1 1 2 2 3 3 4 4 5 5 6 6]),'-','_')       ,...
                     {'HbOrelative', 'HbRrelative', 'HbTrelative', 'HbOrelative', 'HbRrelative', 'HbTrelative'},...
                     {'HbOrelative', 'HbRrelative', 'HbTrelative', 'HbOrelative', 'HbRrelative', 'HbTrelative', 'HbOrelative', 'HbRrelative', 'HbTrelative'}};
plotIdx           = {1:5         ,...
                    repmat([intersect(ExperimentIdx,[2:4,6:8]);intersect(ExperimentIdx,[1 5])],6,1),... 
                    [repmat([1,2],3,1); repmat([3,4],3,1)] ,...
                    [repmat([3 3 5],3,1);repmat([2 4 6],3,1);repmat([1 1 1],3,1)] };
                                    
if isempty(plotIdx{2});plotIdx{2} = repmat([1:8],12,1);end %if the experiment idex does not contain any of the Vazquez experiments plot data for all vazquez experiments. 
%%
plotStructureFieldNames = fieldnames(PlotStructure);
figureHandles = [BOLD_fMRI_figure, Vazquez_HbX_figure, Moon_HbX_figure,Lee_HbX_figure ];


%% plot the data
errorBarHandle = mat2cell(gobjects(sum(cellfun(@numel,plotIdx,'uni',1)),1),cellfun(@numel,plotIdx,'uni',1),[1]);
errorBarHandle = cellfun(@reshape,errorBarHandle,cellfun(@(x)fliplr(size(x)),plotIdx,'uni',0).','uni',0);
for i = 1:numel(figureHandles) % loop over the three main figures
    currentAxes = flipud(findobj(get(figureHandles(i),'Children'), '-depth', 1, 'type', 'axes')); % get all axis objects from figureHandel i 
    for j = 1:numel(currentAxes) % loop over the mumber of axes in figure i
        line(currentAxes(j),currentAxes(j).XLim,[0 0],'LineWidth',0.5,'Color',[0.8 0.8 0.8])% add a thin line along y=0 to indicate baseline
        count_k = 1;
        for k = plotIdx{i}(j,:) % loop over the by plotIdx indicated experiments of PlotStructure field i.
            if ~and(i==2, j>10)
                currentDataStructure = DataStructure.(PlotStructure.(plotStructureFieldNames{i}){k}).(plotVariableNames{i}{j}); % get the data strucute that corresponts to plotStructure: field i, entry k, varaible j. 
                hold(currentAxes(j),'on')
                errorBarHandle{i}(count_k,j) = errorbar(currentAxes(j),currentDataStructure.Time, currentDataStructure.Mean, currentDataStructure.SEM, ' *', 'Color',Colour{i}(k,:),'LineWidth',1, 'CapSize', 2, 'MarkerSize', 4);
                count_k = count_k + 1;
            end
        end
    end
end

%% if simulation Uncertainty should be considered load these results 
plotsimulationUncertainty = 1;

if plotsimulationUncertainty
    load(['.',filesep,'Parameters',filesep,'BestFitsAndParamUC.mat'],'UC_Results')
    if ~isempty(regexp(modelName,'v2[0-9]', 'once'))
        simulationUncertainty = UC_Results.MoonVazquez_Sten2023_v21.simulationUncertainty;
    elseif ~isempty(regexp(modelName,'v3[0-9]', 'once'))
        simulationUncertainty = UC_Results.Lee_Sten2023_v31.simulationUncertainty;
    end
end

%% plot the simulations
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
                if any(variableIdx)% if the logical index has any ture values plot the simualtion i, variable k in plot j 
                    plot(currentAxes(variableIdx),simulation.(simulationFieldNames{i}).t,simulation.(simulationFieldNames{i}).y(:,modelObservablesIndex{i}(k)),'k-','LineWidth',3, 'Color',Colour{j}(simulationFieldInPlotFieldIdx,:));
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
legendStrings ={{'Excitatory 1Hz';
                'Excitatory 20Hz';
                'Forepaw 4Hz'    ;
                'Inhibitory 1Hz' ;
                'Inhibitory 20Hz'};

               {'Excitatory FL'  ;
                'Excitatory 2ms' ;
                'Excitatory 10ms';
                'Excitatory 30ms';
                'Inhibitory FL'  ;
                'Inhibitory 2ms' ;
                'Inhibitory 10ms';
                'Inhibitory 30ms'};

               {'Inhibitory 1Hz' ;
                'Inhibitory 20Hz'};

               {'Whiskers 5Hz ';
                '';
                'SST (SOM) 20Hz ';
                '';
                'nNOS (NO) 20Hz'}};

legendIdex = {[1:5],unique([[1 3 5 8], plotIdx{2}(1,:)],'stable'),[1 2], [1 3 5]};

for i = 1:numel(figureHandles)
    currentAxes = findobj(get(figureHandles(i),'Children'), '-depth', 1, 'type', 'axes');
    count_k = 1;
    for k = legendIdex{i}(1,:)
        legLineHandle{i}(count_k) = line(currentAxes(1),[-11 -10.5], [0 0],'LineWidth',5,'Color',[Colour{i}(k,:)]);
        count_k = count_k + 1;
    end
    legend(legLineHandle{i},legendStrings{i}(legendIdex{i}(1,:)),'AutoUpdate','off','FontSize',globalFontSize)
end

%% set the y axis tick labels for the insert axes 
for i = 1:6
    if i<6; Vazquez_HbX_insertAx(i).YTickLabel = [{num2str(Vazquez_HbX_insertAx(i).YTick(1))} ,repmat({''},1,numel(Vazquez_HbX_insertAx(i).YTick)-2), {num2str(Vazquez_HbX_insertAx(i).YTick(end))}]; end
    if i <=3
        Lee_HbX_insertAx(i).YTickLabel = [{num2str(Lee_HbX_insertAx(i).YTick(1))} ,repmat({''},1,numel(Lee_HbX_insertAx(i).YTick)-2), {num2str(Lee_HbX_insertAx(i).YTick(end))}];
    end 
end

%% plot the line for the stimulation periond after everything else to align with the x-axis
for i = 1:numel(figureHandles)% loop over the three main figures
    currentAxes = flipud(findobj(get(figureHandles(i),'Children'), '-depth', 1, 'type', 'axes')); % get all axis objects from figureHandel i 
    for j = 1:numel(currentAxes) % loop over the mumber of axes in figure i
        k = plotIdx{i}(j,end);
        line(currentAxes(j),[0 DataStructure.(PlotStructure.(plotStructureFieldNames{i}){k}).Stimulation.Length],repmat(currentAxes(j).YLim(1),1,2),'LineWidth',5,'Color',[0 0 0])% add a thin line at thebottomof the plo to illustrate stimulation time
    end
end

%% add A-F sub figure notations to the figures 
addSubFigureNotations = 1;
figLabelFontsize = 22;

if addSubFigureNotations
    setSubFigureNotations(Moon_HbX_ax, figLabelFontsize);
    setSubFigureNotations(Vazquez_HbX_ax([1:6]), figLabelFontsize);
    setSubFigureNotations(Lee_HbX_ax([1:6]), figLabelFontsize);
end

%% print the Chi-threshold
nDataPoints = 0;
for i = 1:numel(ExperimentIdx)
     for j=1:numel(variableNames{i})
         % For variable j in Experiment i add the number of data points for a time point > 0
        nDataPoints =  nDataPoints + nnz(DataStructure.(experimentFieldNames{ExperimentIdx(i)}).(variableNames{i}{j}).Time>0);
     end
end

DoF = nDataPoints;
fprintf('Chi2-threshold for %d degrees of Freedom:\t%0.4f\n',DoF, chi2inv(0.95,DoF))
%% list the parameter values
parameterNames = model.parameters; % put the parameter names into a seperate vector 

% configure the total parameter estimation vector to contian differnt parameter sets for differnt experiments 
for i = 1:size(parameterMappingTable,1) %loop throught the parameter Mapping table.
    parameterNames(parameterMappingTable.TotalParameterIndex{i}) = parameterNames(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.

end

maxNParametersInSetup = max(max(cell2mat(parameterMappingTable.TotalParameterIndex.')),numel(model.parameters)); %compare and get the larger numebr of the largest index in the totalParameter Idx or the number of parameters in the model.
nTheta = maxNParametersInSetup;
nParameterNames = numel(parameterNames);

parameterNameIdx = [true(nTheta - nParameterNames,1);false(2.*nParameterNames - nTheta,1)]; % since we dont have any duplicate parameter values yet this idx vector will just poitn to the (nTheta - nParameterNames) first parameters 
if nTheta > nParameterNames % if the number of parameter values are larger than the numebr of parameters in the model file
    parameterNames(end+1:end+(nTheta - nParameterNames)) = strcat(parameterNames(parameterNameIdx),'_2');
end

parameterTable = table([1:nTheta].',parameterNames.',theta.',10.^(theta.'),'VariableNames',{'Parameter Index','Parameter Name','parameter Value log10',' Parameter value'});
disp (parameterTable);

end
