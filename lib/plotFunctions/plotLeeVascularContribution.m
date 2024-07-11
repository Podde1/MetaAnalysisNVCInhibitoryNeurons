function [] = plotLeeVascularContribution()

modelName = 'Sten2023_v31'; 
ExperimentIdx = 16:21;
VariableNames = {'HbO', 'HbR', 'HbT', 'BOLD', 'CBF'};
load("Parameters\Sten2023_v31_Experiment_16_17_18_19_20_21\fval_(1013.46)_231013_155555_Sten2023_v31.mat", "theta")

[circuitConstants, DataStructure, model, variableNames, ~, modelObservablesIndex, parameterMappingTable, options] = modelSetup(modelName,ExperimentIdx,VariableNames); 

maxNParametersInSetup = max(max(cell2mat(parameterMappingTable.TotalParameterIndex.')),numel(model.parameters)); % compare and get the larger numebr of the largest index in the totalParameter Idx or the number of parameters in the model.
OptParamIndexes = 1:maxNParametersInSetup; % set index of the parameters to be included in parameter estimation. Default is the full maximal number of parameters in the current setup , depending on which experiments are choosen. 
OptParamBool(OptParamIndexes) = true;  % set parameters to optimise to have a true flag
thetaConstants = theta(~OptParamBool); % build the parameter vector that is keept constant during parmater estimation.

[~,simulation] =  objectiveFunction(theta, modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure,ExperimentIdx,variableNames,modelObservablesIndex, parameterMappingTable, model, options, 0, []);

SimSOM16s = simulation.Lee2020_SST_20Hz_16s; 

%parameters
pky1 = 10^theta(ismember(model.parameters, 'ky1'));
pky2 = 10^theta(ismember(model.parameters, 'ky2'));
pky3 = 10^theta(ismember(model.parameters, 'ky3'));

pkNO = 10^theta(ismember(model.parameters, 'kNO'));
pKm_NO = 10^theta(ismember(model.parameters, 'Km_NO'));
pkNO2 = 10^theta(ismember(model.parameters, 'kNO2'));
pKm_NO2 = 10^theta(ismember(model.parameters, 'Km_NO2'));

pVmax_NPY = 10^theta(ismember(model.parameters, 'Vmax_NPY'));
pKm_NPY = 10^theta(ismember(model.parameters, 'Km_NPY'));
pVmax_NPY2 = 10^theta(ismember(model.parameters, 'Vmax_NPY2'));
pKm_NPY2 = 10^theta(ismember(model.parameters, 'Km_NPY2'));

%VSM reactions
NO_in = pkNO.*(SimSOM16s.x(:,18)./(pKm_NO + SimSOM16s.x(:,18)));
SOMNO_in = pkNO2.*(SimSOM16s.x(:,16)./(pKm_NO2 + SimSOM16s.x(:,16)));
NPY_in = pVmax_NPY.*SimSOM16s.x(:,20)./(pKm_NPY+SimSOM16s.x(:,20));
SOMNPY_in = pVmax_NPY2.*SimSOM16s.x(:,17)./(pKm_NPY2 + SimSOM16s.x(:,17));

kvotNO = SOMNO_in./(NO_in+SOMNO_in);
kvotNPY = SOMNPY_in./(NPY_in+SOMNPY_in);

NO_VSM = (1-kvotNO).*SimSOM16s.x(:,19);
SOMNO_VSM = kvotNO.*SimSOM16s.x(:,19);
NPY_VSM = (1-kvotNPY).*SimSOM16s.x(:,21);
SOMNPY_VSM = kvotNPY.*SimSOM16s.x(:,21);

NOvaso = pky1.*(NO_VSM - NO_VSM(1));
SOMNOvaso = pky1.*( SOMNO_VSM - SOMNO_VSM(1));
NPYvaso = pky3.*( NPY_VSM - NPY_VSM(1));
SOMNPYvaso = pky3.*( SOMNPY_VSM - SOMNPY_VSM(1));
PGE2vaso = pky2.*(SimSOM16s.x(:,15)-SimSOM16s.x(1,15));
totalVSM = SimSOM16s.y(:,9);

vasocontributions.PGE2vaso = PGE2vaso;
vasocontributions.NOvaso = NOvaso;
vasocontributions.NPYvaso = -NPYvaso;
vasocontributions.SOMNOvaso = SOMNOvaso;
vasocontributions.SOMNPYvaso = -SOMNPYvaso;
vasocontributions.totalVSM = totalVSM;
vasoNames = fieldnames(vasocontributions);
vasoColor = [[0 0 1]; [0 1 0]; [1 0 0]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]; [0 0 0]];

%% plot structured figure activity to vasular effect
N_state=[];
N_state.N_Pyr = (SimSOM16s.x(:,3) - SimSOM16s.x(1,3))./max(abs(SimSOM16s.x(:,3) - SimSOM16s.x(1,3)));
N_state.N_NO  = (SimSOM16s.x(:,1) - SimSOM16s.x(1,1))./max(abs(SimSOM16s.x(:,1) - SimSOM16s.x(1,1)));
N_state.N_NPY = (SimSOM16s.x(:,2) - SimSOM16s.x(1,2))./max(abs(SimSOM16s.x(:,2) - SimSOM16s.x(1,2)));
N_state.N_SOM = (SimSOM16s.x(:,5) - SimSOM16s.x(1,5))./max(abs(SimSOM16s.x(:,5) - SimSOM16s.x(1,5)));
NNames = fieldnames(N_state);
NColor = [[84 84 255]./255; [121 247 121]./255; [247 77 77]./255; [0.5 0.5 0.5]];

vaso_ymin = min([vasocontributions.NOvaso; vasocontributions.NPYvaso; vasocontributions.PGE2vaso; vasocontributions.SOMNOvaso; vasocontributions.SOMNPYvaso]);
vaso_ymax = max([vasocontributions.NOvaso; vasocontributions.NPYvaso; vasocontributions.PGE2vaso; vasocontributions.SOMNOvaso; vasocontributions.SOMNPYvaso]);

figure('Name','NeuralToVaso', 'WindowState','maximize')
for i = 1:2
    for j = 1:4
        subplot(3, 4, (4*(i-1) + j))
        hold on
        if isequal(i,1)
            plot(SimSOM16s.t, N_state.(NNames{j}), 'LineWidth', 2, 'Color', NColor(j,:))
            ylabel('Normalized neural activity')
            xlim([0 SimSOM16s.t(end)])
            legend(NNames{j}, 'Interpreter','none')
        else
            if isequal(j,4)
                subplot(3, 8, 2*(4*(i-1) + j)-1)
                plot(SimSOM16s.t, vasocontributions.(vasoNames{j}), 'LineWidth', 2, 'Color', vasoColor(j,:))
                ylim([vaso_ymin vaso_ymax])
                xlim([0 SimSOM16s.t(end)])
                ylabel('vasoactive effect')
                xlabel('time (s)')
                legend(vasoNames{j}, 'Interpreter','none')

                subplot(3, 8, 2*(4*(i-1) + j))
                plot(SimSOM16s.t, vasocontributions.(vasoNames{j+1}), 'LineWidth', 2, 'Color', vasoColor(j+1,:))
                legend(vasoNames{j+1}, 'Interpreter','none')
            else
                plot(SimSOM16s.t, vasocontributions.(vasoNames{j}), 'LineWidth', 2, 'Color', vasoColor(j,:))
                legend(vasoNames{j}, 'Interpreter','none')
            end
            ylim([vaso_ymin vaso_ymax])
            xlim([0 SimSOM16s.t(end)])
            ylabel('vasoactive effect')
        end
        xlabel('time (s)')
    end
end
%plot partial VSM contributions 
subplot(3, 4, 9)
plot(SimSOM16s.t, vasocontributions.SOMNOvaso +vasocontributions.SOMNPYvaso, 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', '--')
xlim([0 SimSOM16s.t(end)])
ylim([-0.05 0.1])
xlabel('time (s)')
ylabel('vasoactive effect')
title('SOM total')

subplot(3, 4, 10)
plot(SimSOM16s.t, vasocontributions.SOMNOvaso +vasocontributions.SOMNPYvaso + vasocontributions.PGE2vaso, 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
hold on
plot(SimSOM16s.t, vasocontributions.PGE2vaso, 'LineWidth', 2,  'Color', [0.75 0.75 0.75], 'LineStyle', '--')
xlim([0 SimSOM16s.t(end)])
ylim([-0.05 0.1])
xlabel('time (s)')
ylabel('vasoactive effect')
title('SOM total + PGE2')

subplot(3, 4, 11)
plot(SimSOM16s.t, vasocontributions.SOMNOvaso +vasocontributions.SOMNPYvaso + vasocontributions.PGE2vaso + vasocontributions.NPYvaso, 'LineWidth', 2,  'Color', [0.25 0.25 0.25])
hold on
plot(SimSOM16s.t, vasocontributions.NPYvaso, 'LineWidth', 2,  'Color', [0.75 0.75 0.75], 'LineStyle', '--')
xlim([0 SimSOM16s.t(end)])
ylim([-0.05 0.1])
xlabel('time (s)')
ylabel('vasoactive effect')
title('SOM total + PGE2 + NPY')

subplot(3, 4, 12)
plot(SimSOM16s.t, vasocontributions.(vasoNames{6}), 'LineWidth', 2, 'Color', vasoColor(6,:))
hold on 
plot(SimSOM16s.t, vasocontributions.NOvaso, 'LineWidth', 2,  'Color', [0.75 0.75 0.75], 'LineStyle', '--')
xlabel('time (s)')
ylabel('vasoactive effect')
xlim([0 SimSOM16s.t(end)])
ylim([-0.05 0.1])
legend(vasoNames{6}, 'Interpreter','none')
title('total')


end
