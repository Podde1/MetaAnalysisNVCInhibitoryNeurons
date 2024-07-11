function [] = calculateAndPlotFigure6AUC()
modelName = 'Sten2023_v21';

load(sprintf('./Parameters/%s/AllSimulations_%s.mat',['simUC_',modelName],modelName),'AUC','GABABaseline','N_NOPeakValue')
load(['.',filesep,'Parameters',filesep,'BestFitsAndParamUC.mat'], 'theta_v21')

%% Set experiment index
ExperimentIdx = [1 4 5 8:15]; % specify the experiment index for whihc experiment to run parameter estimation for.
variables2Fit = {'CBF','BOLD','HbT','HbO','HbR'};

%% run model setup
[circuitConstants, DataStructure, modelPropertyNames, variableNames, ~, ~, parameterMappingTable, options] = modelSetup(modelName,ExperimentIdx,variables2Fit);

%% Set up time vectors
stimulationStartTime = 0; % set the stimulation start time. Used to ignore any data points prior to stimulation
experimentNames = fieldnames(DataStructure);% get the field names from the data structure
ExperimentCount = 1;
theta = theta_v21;

for i = ExperimentIdx
    stimulationInfo = DataStructure.(experimentNames{i}).Stimulation; % Used to split the simulation into stimuli and post stimuli simulaiton
    %% set up the time vector for simulating experimetn i
    simulationTime=[];
    for j = 1:numel(variableNames{ExperimentCount})  % loop over the variables to fit in Experiment i.
        simulationTime = [simulationTime,DataStructure.(experimentNames{i}).(variableNames{ExperimentCount}{j}).Time]; % gather all the different time vectors into one
    end

    simulationTime = unique([simulationTime(simulationTime>stimulationStartTime),0, [simulationTime(find(simulationTime>0,1,'first')):0.001:stimulationInfo.Length]]);

    %% set  up the parameter vector for simulating experiment i
    % calculate the number of parameters in the model based on the entries of the parameterMapping table. i.e. the number of model parameters is excatly 1 lower thatn the smallset value in the TotalParameterIndex-column that is also not in the ModelParameterIndex-column.
    nModelParameters = min(setdiff(cell2mat(parameterMappingTable.TotalParameterIndex.'),cell2mat(parameterMappingTable.ModelParameterIndex.')))-1;
    if isempty(nModelParameters); nModelParameters = numel(theta);end % if there is no set difference between ModelParameterIndex and the TotalparameterIndex, get the length of the model parameter vector from theta.
    modelTheta = theta(1:nModelParameters);

    % assign the parameters to use for the model simulation based on the row of the parameterMappingTable that maches the i:tn experiment name.
    modelTheta(parameterMappingTable(experimentNames{i},:).ModelParameterIndex{1}) = theta(parameterMappingTable(experimentNames{i},:).TotalParameterIndex{1});

    % run a single simulation for Experiment i
    simulation.(experimentNames{i}) = simulateSquareStenNVCmodel(modelTheta,simulationTime, modelName, stimulationInfo, circuitConstants, DataStructure.(experimentNames{i}).input, modelPropertyNames,options);

    ExperimentCount = ExperimentCount + 1;
end

%% plot
experiment = {'Moon2021_Inhibitory_1Hz_5s','Moon2021_Inhibitory_1Hz_20s';
    'Moon2021_Inhibitory_20Hz_5s','Moon2021_Inhibitory_20Hz_20s'};
variableIdx = [5,13,15,1];
N = numel(variableIdx) + 2;
Colour = [0.8500  0.3250  0.0980;
          0.9290, 0.6940, 0.1250];

%%
f2 =figure('Units','normalized','WindowState','maximize');
for j = 1:numel(variableIdx)
    for i = 1:size(experiment,2)
        Y = [mean(AUC.(experiment{1,i}).x(:,variableIdx(j))),mean(AUC.(experiment{2,i}).x(:,variableIdx(j)))];

        ax(j,i) = subplot(N,4,i+(j-1)*4);
        hold(ax(j,i),'on')
        for k = 1:2
            barHandle(j,((2*k)-1)+(i-1)) = bar([k],Y(k));
            errorbar(ax(j,i),k,Y(k),...
                [Y(k) - min(AUC.(experiment{k,i}).x(:,variableIdx(j)))],...
                [max(AUC.(experiment{k,i}).x(:,variableIdx(j))) - Y(k)],...
                '.k','LineWidth',2)
        end

        ax(j,i+2) = subplot(N,4,(i+2)+(j-1)*4);
        hold(ax(j,i+2),'on')
        histogram(AUC.(experiment{1,i}).x(:,variableIdx(j)),'FaceColor',Colour(1,:))
        histogram(AUC.(experiment{2,i}).x(:,variableIdx(j)),'FaceColor',Colour(2,:))
        legend(strrep({experiment{1,i},experiment{2,i}},'_',' '))
    end
    set(ax(j,1).YLabel,{'String'},{[strrep(modelPropertyNames.states(variableIdx(j)),'_',' '),' AUC']})
end
%% N_NO peak value
j=5;
for i = 1:size(experiment,2)
    Y = mean(N_NOPeakValue(:,[(2*i)-1 2*i]),1);

    ax(j,i) = subplot(N,4,i+(j-1)*4);
    hold(ax(j,i),'on')
    for k = 1:2
        barHandle(j,((2*k)-1)+(i-1)) =  bar([k],Y(k));

        errorbar(ax(j,i),k,Y(k),...
            Y(k) - min(N_NOPeakValue(:,((2*i)-1)+(k-1)),[],1),...
            max(N_NOPeakValue(:,((2*i)-1)+(k-1)),[],1) - Y(k),...
            '.k','LineWidth',2)
    end
    ax(j,i+2) = subplot(N,4,(i+2)+(j-1)*4);
    hold(ax(j,i+2),'on')
    histogram(N_NOPeakValue(:,((2*i)-1)),'FaceColor',Colour(1,:))
    histogram(N_NOPeakValue(:,((2*i))),'FaceColor',Colour(2,:))
    legend(strrep({experiment{1,i},experiment{2,i}},'_',' '))


end
set(ax(j,1).YLabel,{'String'},{[strrep(modelPropertyNames.states(1),'_',' '),' Peak value']})
%% GABA baseline
j=6;
for i = 1:size(experiment,2)
    Y = mean(GABABaseline(:,[(2*i)-1 2*i]),1);

    ax(j,i) = subplot(N,4,i+(j-1)*4);
    hold(ax(j,i),'on')
    for k = 1:2
        barHandle(j,((2*k)-1)+(i-1)) =  bar([k],Y(k));

        errorbar(ax(j,i),k,Y(k),...
            Y(k) - min(GABABaseline(:,((2*i)-1)+(k-1)),[],1),...
            max(GABABaseline(:,((2*i)-1)+(k-1)),[],1) - Y(k),...
            '.k','LineWidth',2)
    end
    ax(j,i+2) = subplot(N,4,(i+2)+(j-1)*4);
    hold(ax(j,i+2),'on')
    histogram(GABABaseline(:,((2*i)-1)),'FaceColor',Colour(1,:))
    histogram(GABABaseline(:,((2*i))),'FaceColor',Colour(2,:))
    legend(strrep({experiment{1,i},experiment{2,i}},'_',' '))
end
set(ax(j,1).YLabel,{'String'},{[strrep(modelPropertyNames.states(5),'_',' '),' Stimulation Baseline']})

set(ax(1,1).Title,{'String'},{'5 sec stimulation'})
set(ax(1,2).Title,{'String'},{'20 sec stimulation'})
set(ax(1,3).Title,{'String'},{'5 sec stimulation'})
set(ax(1,4).Title,{'String'},{'20 sec stimulation'})

set(barHandle(:,[1,2]),{'FaceColor'},{Colour(1,:)})
set(barHandle(:,[3,4]),{'FaceColor'},{Colour(2,:)})

set(ax,{'FontSize'},{18})
set(ax(:,[1 2]),{'XTick','XTickLabel'},{[1 2],{'1 Hz','20 Hz'}})

f2.Children([3 5 9 11 15 17 21 23 27 29 33 35]) = f2.Children([5 3 11 9 17 15 23 21 29 27 35 33]); % Restructure the order of the f2. Children such that the maximuise function does no reorder the axes.
pause(3); % a pause is needed for Matlab to render th axes before resizing them 
maximizeAxes(f2,6,4);
setSubFigureNotations(reshape(ax.',24,1),24);

%%
variableIdx = [5 1];
figure('Units','normalized','WindowState','maximize')
for i = 1:2
    
    for j =[2 1]
    ax2(((2*i)-1)+(j-1)) = subplot(4,1,((2*i)-1)+(j-1));
    hold on
    plot(simulation.(experiment{j,2}).t,simulation.(experiment{j,2}).x(:,variableIdx(i))-simulation.(experiment{j,2}).x(1,variableIdx(i)), 'LineWidth',2,'Color',Colour(j,:))
    line([-1 70],[0 0],'Linewidth',1,'Color',[0.5 0.5 0.5]);

    set(ax2(((2*i)-1)+(j-1)).YLabel,{'String'},{strrep(modelPropertyNames.states(variableIdx(i)),'_',' ')})
    set(ax2(((2*i)-1)+(j-1)).XLabel,{'String'},{'Time (s)'})
    end
end
set(ax2,{'FontSize','XLim'},{18,[-1 30]})

%% 
Vaz_experiment = {'Vazques2018_Excitatory_30ms','Vazques2018_Inhibitory_30ms'};
Colour2 = [0       0.2627  0.6510 ;
           0.8500  0.3250  0.0980];
f3 =figure('Units','normalized','WindowState','maximize');   
    Y = [mean(AUC.(Vaz_experiment{1}).x(:,[13,15]));mean(AUC.(Vaz_experiment{2}).x(:,[13,15]))];
    Y_diff = [Y(1,:) - min(AUC.(Vaz_experiment{1}).x(:,[13,15])),Y(2,:) - min(AUC.(Vaz_experiment{2}).x(:,[13,15]));
              max(AUC.(Vaz_experiment{1}).x(:,[13,15])) - Y(1,:),max(AUC.(Vaz_experiment{2}).x(:,[13,15])) - Y(2,:)];

    barHandle2 = bar([1 2],Y);
    hold on

for i =1:numel(Vaz_experiment)
    errorbar(barHandle2(i).XEndPoints,Y(:,i).',...
            Y_diff(1,[i i+2]),...
            Y_diff(2,[i i+2]),...
            '.k','LineWidth',2)
    set(barHandle2(i),{'FaceColor'},{Colour2(i,:)});
end
ax3=gca;
set(ax3,{'FontSize','YLim','XTickLabel'},{18,[0 0.35],{'Excitatory','Inhibitory'}})
ylabel('AUC')
legend(barHandle2,{'PGE_{2,VSM}','NO_{VSM}'})


%% vsm simulations
pky1 = 10^theta(ismember(modelPropertyNames.parameters, 'ky1'));
pky2 = 10^theta(ismember(modelPropertyNames.parameters, 'ky2'));

vasoColor = [[0 1 0]; [0 0 1]; [0.5 0.5 0.5]];

figure()
dpi = 96;
a4WidthPixels = round(0.4*210 / 25.4 * dpi);
a4HeightPixels = round(0.4*297 / 25.4 * dpi);
set(gcf, 'Position', [0 0 a4WidthPixels a4HeightPixels]);
sgtitle('VSM effect')

fnames = experiment(3:4)';
for i = 1:length(fnames)
    NOvsm   =  pky1.*(simulation.(fnames{i}).x(:,15) - simulation.(fnames{i}).x(1,15));
    PGE2vsm =  pky2.*(simulation.(fnames{i}).x(:,13) - simulation.(fnames{i}).x(1,13));
    total   =  NOvsm + PGE2vsm;

    nexttile
    hold on
    plot(simulation.(fnames{i}).t([1 end]), [0 0],   'LineWidth', 1, 'Color', [0 0 0])
    plot(simulation.(fnames{i}).t, NOvsm,   'LineWidth', 2, 'Color', vasoColor(1,:))
    plot(simulation.(fnames{i}).t, PGE2vsm, 'LineWidth', 2, 'Color', vasoColor(2,:))
    plot(simulation.(fnames{i}).t, total,   'LineWidth', 1, 'Color', vasoColor(3,:))
    ylabel('vsm effect (AU)')
    xlabel('time (s)')
    xlim([0 25])
    ylim([-1.5 1.5])

    ax = gca;
    set(ax.XAxis, 'FontSize', 10)
    set(ax.YAxis, 'FontSize', 10)
end

end