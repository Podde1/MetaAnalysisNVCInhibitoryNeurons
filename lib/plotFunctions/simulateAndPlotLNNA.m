function [] = simulateAndPlotLNNA()

modelName = 'Sten2023_v21';

%% run model setup 
ExperimentIdx = [5,8];
variables2Fit = {'CBF','BOLD','HbT','HbO','HbR'};

[circuitConstants, DataStructure, modelPropertyNames, ~, ~, ~, ~, options] = modelSetup(modelName,ExperimentIdx,variables2Fit);

%% load NOS blocker data
load(['.',filesep,'Data',filesep,'Vazquez_CerebCortex_2018_Data.mat'],'DataVaquez2018')
%%
simulationTime = [DataVaquez2018.Experiment4_L_NNA_30ms.CBF.Time, DataVaquez2018.Experiment4_L_NNA_FL.CBF.Time];
stimulationStartTime = 0;
simulationTime = unique([simulationTime(simulationTime>stimulationStartTime),0]); % sort out duplicate values and sort the vector. Also sort out values <0. add 0 such that this time poind can be includend in plots

stimulationInfo_inhbi = DataStructure.Vazques2018_Inhibitory_30ms.Stimulation; % Used to split the simulation into stimuli and post stimuli simulaiton
stimulationInfo_sensory = DataStructure.Vazques2018_Inhibitory_FL.Stimulation;

inputTable_inhib = DataStructure.Vazques2018_Inhibitory_30ms.input;
inputTable_Sensory = inputTable_inhib;
inputTable_Sensory.Stimulation(3)=1;
%% load parametervalues 
load(['.',filesep,'Parameters',filesep,'BestFitsAndParamUC.mat'],'UC_Results', 'theta_v21')

theta_v21(find(strcmp(modelPropertyNames.parameters,'kNOS'))) = -8; % set kNOS parameter to minimum value

modelTheta_inhib_30ms = theta_v21(1:numel(modelPropertyNames.parameters));
modelTheta_sensory = theta_v21(1:numel(modelPropertyNames.parameters));

modelTheta_inhib_30ms(UC_Results.MoonVazquez_Sten2023_v21.parameterMappingTable(4,:).ModelParameterIndex{1}) = theta_v21(UC_Results.MoonVazquez_Sten2023_v21.parameterMappingTable(4,:).TotalParameterIndex{1});
modelTheta_sensory(UC_Results.MoonVazquez_Sten2023_v21.parameterMappingTable(3,:).ModelParameterIndex{1}) = theta_v21(UC_Results.MoonVazquez_Sten2023_v21.parameterMappingTable(1,:).TotalParameterIndex{1});

%%
simulation_Inhib_30ms = simulateSquareStenNVCmodel(modelTheta_inhib_30ms.',simulationTime, modelName, stimulationInfo_inhbi, circuitConstants,  inputTable_inhib, modelPropertyNames,options);
simulation_Sensory = simulateSquareStenNVCmodel(modelTheta_sensory.',simulationTime, modelName, stimulationInfo_sensory, circuitConstants,  inputTable_Sensory, modelPropertyNames,options);

load(['.',filesep,'Parameters',filesep,'LNNA_simUC.mat'], "CBF_LNNA_sensory","CBF_LNNA_inhib")

sensorysimulation = [CBF_LNNA_sensory.min;CBF_LNNA_sensory.max;simulation_Sensory.y(:,1).'];
inhibitorysimulation = [CBF_LNNA_inhib.min;CBF_LNNA_inhib.max;simulation_Inhib_30ms.y(:,1).'];

AUC_Sensory_data = trapz(DataVaquez2018.Experiment4_L_NNA_FL.CBF.Mean + [-1;1;0].*DataVaquez2018.Experiment4_L_NNA_FL.CBF.SEM,2);
AUC_inhib_data = trapz(DataVaquez2018.Experiment4_L_NNA_10ms.CBF.Mean + [-1;1;0].*DataVaquez2018.Experiment4_L_NNA_10ms.CBF.SEM,2);

AUC_inhib_data_noBlock = trapz(DataVaquez2018.Experiment5_Inhibitory_30ms.CBF.Mean + [-1; 1; 0].* DataVaquez2018.Experiment5_Inhibitory_30ms.CBF.SEM,2);
AUC_Sensory_data_noBlock = trapz(DataVaquez2018.Experiment5_Inhibitory_FL.CBF.Mean + [-1; 1; 0].* DataVaquez2018.Experiment5_Inhibitory_FL.CBF.SEM,2);

AUC_sensory_sim = trapz(sensorysimulation,2);
AUC_inhib_sim = trapz(inhibitorysimulation,2);

peakValue_inhib_simulation = max(inhibitorysimulation,[],2);
peakValue_sensory_simulation = max(sensorysimulation,[],2);

peakValue_Sensory_data = max(DataVaquez2018.Experiment4_L_NNA_FL.CBF.Mean + [-1;1;0].*DataVaquez2018.Experiment4_L_NNA_FL.CBF.SEM,[],2);
peakValue_inhib_data = max(DataVaquez2018.Experiment4_L_NNA_10ms.CBF.Mean + [-1;1;0].*DataVaquez2018.Experiment4_L_NNA_10ms.CBF.SEM,[],2);

peakValue_inhib_data_noBlock = max(DataVaquez2018.Experiment5_Inhibitory_30ms.CBF.Mean + [-1; 1; 0].* DataVaquez2018.Experiment5_Inhibitory_30ms.CBF.SEM,[],2);
peakValue_Sensory_data_noBlock = max(DataVaquez2018.Experiment5_Inhibitory_FL.CBF.Mean + [-1; 1; 0].* DataVaquez2018.Experiment5_Inhibitory_FL.CBF.SEM,[],2);

%% plot AUC values
figure()
ax2(1) = subplot(1,2,1);
bar([1:3],[AUC_Sensory_data_noBlock(3),AUC_Sensory_data(3),AUC_sensory_sim(3)],'FaceColor',[0.5  0.5  0.5])
title('Sensory Stimulation')
set(ax2(1),{'XTickLabel'},{{'Standard data','NO blocked data','NO Blocked simulation'}})

ax2(2) = subplot(1,2,2);
bar([1:3],[AUC_inhib_data_noBlock(3),AUC_inhib_data(3),AUC_inhib_sim(3)],'FaceColor',[0.8500  0.3250  0.0980])
title('Inhibitory Stimulation')
set(ax2(2),{'XTickLabel'},{{'Standard data','NO blocked data','NO Blocked simulation'}})

set(ax2,{'Ylim'},{[0,450]});
sgtitle('AUC values')

%% plot peak values
figure()
ax3(1) = subplot(1,2,1);
hold on
bar([1:3],[peakValue_Sensory_data_noBlock(3),peakValue_Sensory_data(3),peakValue_sensory_simulation(3)],'FaceColor',[0.5  0.5  0.5])
errorbar([1:3],[peakValue_Sensory_data_noBlock(3),peakValue_Sensory_data(3),peakValue_sensory_simulation(3)],[peakValue_Sensory_data_noBlock(3),peakValue_Sensory_data(3),peakValue_sensory_simulation(3)]-[peakValue_Sensory_data_noBlock(1),peakValue_Sensory_data(1),peakValue_sensory_simulation(1)],[peakValue_Sensory_data_noBlock(2),peakValue_Sensory_data(2),peakValue_sensory_simulation(2)]-[peakValue_Sensory_data_noBlock(3),peakValue_Sensory_data(3),peakValue_sensory_simulation(3)],' .k');
title('Sensory Stimulation')
set(ax3(1),{'XTick','XTickLabel'},{[1:3],{'Standard data','NO blocked data','NO Blocked simulation'}})

ax3(2) = subplot(1,2,2);
hold on
bar([1:3],[peakValue_inhib_data_noBlock(3),peakValue_inhib_data(3),peakValue_inhib_simulation(3)],'FaceColor',[0.8500  0.3250  0.0980])
errorbar([1:3],[peakValue_inhib_data_noBlock(3),peakValue_inhib_data(3),peakValue_inhib_simulation(3)],[peakValue_inhib_data_noBlock(3),peakValue_inhib_data(3),peakValue_inhib_simulation(3)]-[peakValue_inhib_data_noBlock(1),peakValue_inhib_data(1),peakValue_inhib_simulation(1)],[peakValue_inhib_data_noBlock(2),peakValue_inhib_data(2),peakValue_inhib_simulation(2)]-[peakValue_inhib_data_noBlock(3),peakValue_inhib_data(3),peakValue_inhib_simulation(3)],' .k')
title('Inhibitoru Stimulation')
set(ax3(2),{'XTick','XTickLabel'},{[1:3],{'Standard data','NO blocked data','NO Blocked simulation'}})

set(ax3,{'Ylim'},{[0,25]});
sgtitle('Peak values')

%% plot simualtion 
figHandle =  figure();

ax1 = subplot(2,1,1);
hold on
errorbar(DataVaquez2018.Experiment4_L_NNA_10ms.CBF.Time,DataVaquez2018.Experiment4_L_NNA_10ms.CBF.Mean,DataVaquez2018.Experiment4_L_NNA_10ms.CBF.SEM,' *k','Color',[0.8500  0.3250  0.0980],'LineWidth',2)
plot(simulation_Inhib_30ms.t,simulation_Inhib_30ms.y(:,1),'k-','LineWidth',3, 'Color',[0.8500  0.3250  0.0980])

fill([simulation_Inhib_30ms.t.',fliplr(simulation_Inhib_30ms.t.')],[CBF_LNNA_inhib.min,fliplr(CBF_LNNA_inhib.max)],'',...
                            'FaceColor',[0.8500  0.3250  0.0980],...
                            'FaceAlpha',0.2, ...
                            'EdgeColor',[0.8500  0.3250  0.0980])
ylabel('\Delta CBF(%)')
xlabel('Time (s)')
line([0 4],[-4 -4],'LineWidth',5,'Color',[0 0 0])
set(ax1,{'FontSize','YLim','XLim'},{24,[-4 15],[-5 55]})
lin(1)= line([-10 -11], [0 0],'LineWidth',5,'Color',[0.8500  0.3250  0.0980]);
lin(2)= line([-10 -11], [0 0],'LineWidth',5,'Color',[0.5     0.5     0.5]);
legend(lin,{'Inhibitory OG stimulation with NOS blocked' ,'Sensory stimulation with NOS blocked'})


ax2 = subplot(2,1,2);
hold on
errorbar(DataVaquez2018.Experiment4_L_NNA_FL.CBF.Time,DataVaquez2018.Experiment4_L_NNA_FL.CBF.Mean,DataVaquez2018.Experiment4_L_NNA_FL.CBF.SEM,' *k','Color',[0.5     0.5     0.5],'LineWidth',2)
fill([simulation_Sensory.t.',fliplr(simulation_Sensory.t.')],[CBF_LNNA_sensory.min,fliplr(CBF_LNNA_sensory.max)],'',...
                            'FaceColor',[0.5     0.5     0.5],...
                            'FaceAlpha',0.2, ...
                            'EdgeColor',[0.5     0.5     0.5])
line([0 4],[-4 -4],'LineWidth',5,'Color',[0 0 0])
set(ax2,{'FontSize','YLim','XLim'},{24,[-4 15],[-5 55]})
ylabel('\Delta CBF(%)')
xlabel('Time (s)')

end