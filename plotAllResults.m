function [] = plotAllResults()
%% This script runs all relevant plot functions to plot all results presented in the paper entitled "". 
% Please ensure that a MATLAB release of R2023a or newer is used to run
% this script 

%Please ensure that the AMICI toolbox is installed  before runnig this script (https://github.com/AMICI-dev/AMICI) 

% Finally, please ensure that all relevant models have been compiled before running this script.  To do this please run the "GenerateModels()" script has been run in
% MATLAB R2017b.
%% set up
close all

addpath(genpath(fullfile(pwd,'./Data')))
addpath(genpath(fullfile(pwd,'./lib')))
addpath(genpath(fullfile(pwd,'./Models')))
addpath(genpath(fullfile(pwd,'./Parameters')))

%% check Compatibility
 if ~checkCompatibility() % if base requirements are not met terminate script. 
     return
 end

%% neuronal contributions - Figures 2, S1, S2, and S3
fprintf('Plotting BOLD contributions \n')
plotBOLDcontributionFromFile()
fprintf('Plotting HbR contributions \n')
plotHbRcontributionFromFile()
fprintf('Plotting CMRO2 contributions \n')
plotCMRO2contributionFromFile()
fprintf('Plotting CBV contributions \n')
plotCBVcontributionFromFile()

%% model estimation - Figures 3, 5, and 6
load(['.',filesep,'Parameters',filesep,'BestFitsAndParamUC.mat'], 'theta_v21')
fprintf('Plotting Vazquez and Moon model estimations \n')
MoonVazFigures = plotCurrenSolution(theta_v21, 'Sten2023_v21', [1 4 5 8:15], {'CBF','BOLD','HbT','HbO','HbR'});
close(MoonVazFigures(4)); % Lee data

% Figure 7
load(['.',filesep,'Parameters',filesep,'BestFitsAndParamUC.mat'], 'theta_v31')
fprintf('Plotting Lee model estimations \n')
LeeFigures = plotCurrenSolution(theta_v31, 'Sten2023_v31', [16:21], {'CBF','BOLD','HbT','HbO','HbR'});
close(LeeFigures(1:3)); % Vazquex and Moon data

%% Qualitative analysis
% Figure 6 analysis 
fprintf('Plotting Moon analysis \n')
calculateAndPlotFigure6AUC()

% Figure S4 analysis
fprintf('Plotting Lee analysis \n')
plotLeeVascularContribution() 

%% model validations 
% NOS experiment Vazquez 2018 - Figure 3G
fprintf('Plotting NOS validation\n')
simulateAndPlotLNNA()

% 2ms and 10ms pulse widhts Vazquez 2018 - Figure 4
fprintf('Plotting Vazquez 2ms and 10ms validation \n')
plotVazquezValidation()

fprintf('\n \n #---------- ALL RESULTS HAVE BEEN PLOTTED ----------# \n')
end