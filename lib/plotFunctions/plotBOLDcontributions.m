close all
clear all 

%% Set up
shouldIPlot = 0; 
modelName = 'Sten2023_v21_BOLDcontributions';
ExperimentIdx = [1 4 5 8 9:15];
VariableNames = {'HbO', 'HbR', 'HbT', 'BOLD', 'CBF'};

[circuitConstants, DataStructure, model, variableNames, saveString, modelObservablesIndex, parameterMappingTable, options] = modelSetup(modelName,ExperimentIdx,VariableNames); 

maxNParametersInSetup = max(max(cell2mat(parameterMappingTable.TotalParameterIndex.')),numel(model.parameters)); %compare and get the larger number of the largest index in the totalParameter Idx or the number of parameters in the model.
OptParamIndexes = 1:maxNParametersInSetup; %set index of the parameters to be included in parameter estimation. Default is the full maximal number of parameters in the current setup , depending on which experiments are choosen. 
OptParamBool(OptParamIndexes) = true; %set parameters to optimise to have a true flag
thetaConstants = []; % build the parameter vector that is keept constant during parmater estimation.

%% set up indexes to various model properties
[~, idxCa_NO]  = ismember('Ca_NO' ,   model.states);
[~, idxCa_NPY] = ismember('Ca_NPY',   model.states);
[~, idxCa_Pyr] = ismember('Ca_Pyr',   model.states);
[~, idxCa_Ast] = ismember('Ca_Astro', model.states);
[~, idxNOvsm]  = ismember('NOvsm',    model.states);
[~, idxNPYvsm] = ismember('NPYvsm',   model.states);
[~, idxPyrvsm] = ismember('PGE2vsm',  model.states);
[~, idxV1]     = ismember('V1',       model.states);
[~, idxV2]     = ismember('V2',       model.states);
[~, idxV3]     = ismember('V3',       model.states);
[~, idxnO2_1]  = ismember('nO2_1',    model.states);
[~, idxnO2_2]  = ismember('nO2_2',    model.states);
[~, idxnO2_3]  = ismember('nO2_3',    model.states);

[~, idxCMRO2]  = ismember('CMRO2',   model.observables );
[~, idxjO2_1]  = ismember('jO2_1',   model.observables );
[~, idxjO2_2]  = ismember('jO2_2',   model.observables );
[~, idxjO2_3]  = ismember('jO2_3',   model.observables );
[~, idxCO2_01] = ismember('CO2_01',  model.observables );

idxHbR         = find(strcmp('HbR',  model.observables));
idxHbR = idxHbR(end);
[~, idxBOLD_OIS] = ismember('BOLD_OIS' , model.observables);

%% Extract parameters
load('Parameters\BestFitsAndParamUC.mat', 'UC_Results', 'theta_v21');

allParams = theta_v21;
count = 1;
MCMCresStruct = UC_Results.MoonVazquez_Sten2023_v21.MCMCUncertainty;;
for field = fieldnames(MCMCresStruct)' 
    for subfield = fieldnames(MCMCresStruct.(field{:}))'
        if contains(subfield{:}, 'param', IgnoreCase=true)
            count = count + 1;
            allParams(count, :) = MCMCresStruct.(field{:}).(subfield{:}); 
        end
    end 
end

%% Initialize tables
allExperiments = fieldnames(DataStructure);
Experiments    = allExperiments(ExperimentIdx);
BOLDTable      = cell(numel(Experiments), 5);
CBVTable       = cell(numel(Experiments), 5);
CMRO2Table     = cell(numel(Experiments), 5);
HbRTable       = cell(numel(Experiments), 5);

for i = 1:size(BOLDTable,1)
    BOLDTable{i,1}  = Experiments{i}; % Name
    CBVTable{i,1}   = Experiments{i}; % Name
    CMRO2Table{i,1} = Experiments{i}; % Name
    HbRTable{i,1}   = Experiments{i}; % Name

    for j = 2:5
        BOLDTable{i,j}  = cell(size(allParams,1),1);
        CBVTable{i,j}   = cell(size(allParams,1),1);
        CMRO2Table{i,j} = cell(size(allParams,1),1);
        HbRTable{i,j}   = cell(size(allParams,1),1);
        
        for k = 1:size(allParams,1)
            BOLDTable{i,j}{k}  = NaN;
            CBVTable{i,j}{k}   = NaN;
            CMRO2Table{i,j}{k} = NaN;
            HbRTable{i,j}{k}   = NaN;
        end
    end
end

%% Simulate the different neuronal contributions 
for idParam = 1:size(allParams, 1)

    if isequal(mod(idParam, 10), 0)
        fprintf("Unceratainty eastimation is %.1f %s done \n", 100*(idParam/size(allParams, 1)), '%')
    end

    theta = allParams(idParam, :);
    
    warning('off')
    [~,simulationBOLD] = objectiveFunction(theta, modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure, ExperimentIdx, variableNames, modelObservablesIndex, parameterMappingTable, model, options, 0, []);
    warning('on')

    % define parameter values 
    pky1         = 10^theta(ismember(model.parameters, 'ky1'));
    pky2         = 10^theta(ismember(model.parameters, 'ky2'));
    pky3         = 10^theta(ismember(model.parameters, 'ky3'));
    kCMRO2_NO    = 10^theta(ismember(model.parameters, 'kCMRO2_NO'));
    kCMRO2_NPY   = 10^theta(ismember(model.parameters, 'kCMRO2_NPY'));
    kCMRO2_Pyr   = 10^theta(ismember(model.parameters, 'kCMRO2_Pyr'));
    kCMRO2_Astro = 10^theta(ismember(model.parameters, 'kCMRO2_Astro'));
    k_yBOLD      = 10^theta(ismember(model.parameters, 'k_yBOLD'));
    k_yBOLD2     = 10^theta(ismember(model.parameters, 'k_yBOLD2'));
    
    for i = 1:length(Experiments)
        Experiment = Experiments{i};
        
        if simulationBOLD.(Experiment).status > 0 % check that simulation did not crash
            % vsm stimulation effect
            stim_circNO   =  pky1.*(simulationBOLD.(Experiment).x(:,idxNOvsm)  - simulationBOLD.(Experiment).x(1,idxNOvsm));
            stim_circPyr  =  pky2.*(simulationBOLD.(Experiment).x(:,idxPyrvsm) - simulationBOLD.(Experiment).x(1,idxPyrvsm));
            stim_circNPY  = -pky3.*(simulationBOLD.(Experiment).x(:,idxNPYvsm) - simulationBOLD.(Experiment).x(1,idxNPYvsm));
            stim_circabs = abs(stim_circNO ) + abs(stim_circNPY) + abs(stim_circPyr) + 1e-10;
            
            % fractional contributions to the total vsm effect 
            NOfrac  = abs(stim_circNO )./stim_circabs;
            NPYfrac = abs(stim_circNPY)./stim_circabs;
            Pyrfrac = abs(stim_circPyr)./stim_circabs; 
            
            M1    = lt(stim_circNO, 0);     M2    = lt(stim_circNPY, 0);     M3    = lt(stim_circPyr, 0); % if NOeffect, NPYeffect, and Pyreffect is negative
            M1inv = ge(stim_circNO, 0);     M2inv = ge(stim_circNPY, 0);     M3inv = ge(stim_circPyr, 0); % if NOeffect, NPYeffect, and Pyreffect is positive
            
            % calculate and balance out the missing volume from the negative contributions
            NO_frac  = (stim_circNO ./stim_circabs)  +    ( (M2.*((NOfrac ./(NOfrac  + M3inv.*Pyrfrac)).*(2*NPYfrac)) + M3.*((NOfrac ./(NOfrac  + M2inv.*NPYfrac)).*(2*Pyrfrac)) ) ).*M1inv;
            NPY_frac = (stim_circNPY./stim_circabs)  +    ( (M1.*((NPYfrac./(NPYfrac + M3inv.*Pyrfrac)).*(2*NOfrac )) + M3.*((NPYfrac./(NPYfrac + M1inv.*NOfrac )).*(2*Pyrfrac)) ) ).*M2inv;
            Pyr_frac = (stim_circPyr./stim_circabs)  +    ( (M1.*((Pyrfrac./(Pyrfrac + M2inv.*NPYfrac)).*(2*NOfrac )) + M2.*((Pyrfrac./(Pyrfrac + M1inv.*NOfrac )).*(2*NPYfrac)) ) ).*M3inv;
            
            NO_frac(1)  = 0;
            NPY_frac(1) = 0;
            Pyr_frac(1) = 0;
            
            V1 = simulationBOLD.(Experiment).x(:,idxV1);
            V2 = simulationBOLD.(Experiment).x(:,idxV2);
            V3 = simulationBOLD.(Experiment).x(:,idxV3);
            
            % calculate the respective volume effects of the neurons 
            V1NO  = V1(1) + NO_frac .*(V1-V1(1));      V1NPY = V1(1) + NPY_frac.*(V1-V1(1));      V1Pyr = V1(1) + Pyr_frac.*(V1-V1(1));    
            V2NO  = V2(1) + NO_frac .*(V2-V2(1));      V2NPY = V2(1) + NPY_frac.*(V2-V2(1));      V2Pyr = V2(1) + Pyr_frac.*(V2-V2(1));
            V3NO  = V3(1) + NO_frac .*(V3-V3(1));      V3NPY = V3(1) + NPY_frac.*(V3-V3(1));      V3Pyr = V3(1) + Pyr_frac.*(V3-V3(1));
            
            %% CMRO2 contribution
            CMRO2 = simulationBOLD.(Experiment).y(:,idxCMRO2);
            CBVss = 1;
            
            % active (over basal) CMRO2 contribution 
            CMRO2NO  = CMRO2(1) .* (kCMRO2_NO    .* (simulationBOLD.(Experiment).x(:,idxCa_NO ) - simulationBOLD.(Experiment).x(1,idxCa_NO ) ));
            CMRO2NPY = CMRO2(1) .* (kCMRO2_NPY   .* (simulationBOLD.(Experiment).x(:,idxCa_NPY) - simulationBOLD.(Experiment).x(1,idxCa_NPY) ));
            CMRO2Pyr = CMRO2(1) .* (kCMRO2_Pyr   .* (simulationBOLD.(Experiment).x(:,idxCa_Pyr) - simulationBOLD.(Experiment).x(1,idxCa_Pyr) ));
            CMRO2Ast = CMRO2(1) .* (kCMRO2_Astro .* (simulationBOLD.(Experiment).x(:,idxCa_Ast) - simulationBOLD.(Experiment).x(1,idxCa_Ast) ));
            CMRO2tot = abs(CMRO2NO) + abs(CMRO2NPY) + abs(CMRO2Pyr) + abs(CMRO2Ast) + 1e-10;
            
            % oxygen extracted from the vessels to tissue
            jO2_1 = simulationBOLD.(Experiment).y(:,idxjO2_1);
            jO2_2 = simulationBOLD.(Experiment).y(:,idxjO2_2);
            jO2_3 = simulationBOLD.(Experiment).y(:,idxjO2_3);
            
            % oxygen that should not have been consumed when isolating the neurons i.e. Only NO should not include NPY, Pyr, and Ast
            jO2_1NORet = jO2_1.*((CMRO2NPY + CMRO2Pyr + CMRO2Ast)./CMRO2);     jO2_1NPYRet = jO2_1.*((CMRO2NO + CMRO2Pyr + CMRO2Ast)./CMRO2);     jO2_1PyrRet = jO2_1.*((CMRO2NO + CMRO2NPY + CMRO2Ast)./CMRO2);     jO2_1AstRet = jO2_1.*((CMRO2NO + CMRO2NPY + CMRO2Pyr)./CMRO2);
            jO2_2NORet = jO2_2.*((CMRO2NPY + CMRO2Pyr + CMRO2Ast)./CMRO2);     jO2_2NPYRet = jO2_2.*((CMRO2NO + CMRO2Pyr + CMRO2Ast)./CMRO2);     jO2_2PyrRet = jO2_2.*((CMRO2NO + CMRO2NPY + CMRO2Ast)./CMRO2);     jO2_2AstRet = jO2_2.*((CMRO2NO + CMRO2NPY + CMRO2Pyr)./CMRO2);
            jO2_3NORet = jO2_3.*((CMRO2NPY + CMRO2Pyr + CMRO2Ast)./CMRO2);     jO2_3NPYRet = jO2_3.*((CMRO2NO + CMRO2Pyr + CMRO2Ast)./CMRO2);     jO2_3PyrRet = jO2_3.*((CMRO2NO + CMRO2NPY + CMRO2Ast)./CMRO2);     jO2_3AstRet = jO2_3.*((CMRO2NO + CMRO2NPY + CMRO2Pyr)./CMRO2);
            
            % oxygen amount
            nO2_1 = simulationBOLD.(Experiment).x(:,idxnO2_1);
            nO2_2 = simulationBOLD.(Experiment).x(:,idxnO2_2);
            nO2_3 = simulationBOLD.(Experiment).x(:,idxnO2_3);
            
            % oxygen amount for isolated CMRO2 effect
            nO2_1no = nO2_1 + jO2_1NORet;        nO2_1npy = nO2_1 + jO2_1NPYRet;        nO2_1pyr = nO2_1 + jO2_1PyrRet;        nO2_1ast = nO2_1 + jO2_1AstRet;    
            nO2_2no = nO2_2 + jO2_2NORet;        nO2_2npy = nO2_2 + jO2_2NPYRet;        nO2_2pyr = nO2_2 + jO2_2PyrRet;        nO2_2ast = nO2_2 + jO2_2AstRet;    
            nO2_3no = nO2_3 + jO2_3NORet;        nO2_3npy = nO2_3 + jO2_3NPYRet;        nO2_3pyr = nO2_3 + jO2_3PyrRet;        nO2_3ast = nO2_3 + jO2_3AstRet;    
            
            % estimated oxygen amount from the vascular activity of the neurons
            nO2_1NO = nO2_1no(2) + NO_frac .*(nO2_1no-nO2_1no(2));       nO2_1NPY = nO2_1npy(2) + NPY_frac .*(nO2_1npy-nO2_1npy(2));       nO2_1Pyr = nO2_1pyr(2) + Pyr_frac .*(nO2_1pyr-nO2_1pyr(2));       nO2_1Ast = nO2_1ast(2);
            nO2_2NO = nO2_2no(2) + NO_frac .*(nO2_2no-nO2_2no(2));       nO2_2NPY = nO2_2npy(2) + NPY_frac .*(nO2_2npy-nO2_2npy(2));       nO2_2Pyr = nO2_2pyr(2) + Pyr_frac .*(nO2_2pyr-nO2_2pyr(2));       nO2_2Ast = nO2_2ast(2);
            nO2_3NO = nO2_3no(2) + NO_frac .*(nO2_3no-nO2_3no(2));       nO2_3NPY = nO2_3npy(2) + NPY_frac .*(nO2_3npy-nO2_3npy(2));       nO2_3Pyr = nO2_3pyr(2) + Pyr_frac .*(nO2_3pyr-nO2_3pyr(2));       nO2_3Ast = nO2_3ast(2);
            
            % calculate the oxygen concentrations 
            CO2_max = 9.26;
            CO2_01 = simulationBOLD.(Experiment).y(:,idxCO2_01);
            
            CO2_12NO = nO2_1NO./V1NO;    CO2_12NPY = nO2_1NPY./V1NPY;    CO2_12Pyr = nO2_1Pyr./V1Pyr;    CO2_12Ast = nO2_1Ast./V1(1);
            CO2_23NO = nO2_2NO./V2NO;    CO2_23NPY = nO2_2NPY./V2NPY;    CO2_23Pyr = nO2_2Pyr./V2Pyr;    CO2_23Ast = nO2_2Ast./V2(1);
            CO2_34NO = nO2_3NO./V3NO;    CO2_34NPY = nO2_3NPY./V3NPY;    CO2_34Pyr = nO2_3Pyr./V3Pyr;    CO2_34Ast = nO2_3Ast./V3(1);
            
            % calculate the oxygen saturation, 95% is set as a pysiological cap
            SaO2NO = min(0.95, ((CO2_01   + CO2_12NO)/2)./(CO2_max));     SaO2NPY = min(0.95, ((CO2_01    + CO2_12NPY)/2)./(CO2_max));     SaO2Pyr = min(0.95, ((CO2_01    + CO2_12Pyr)/2)./(CO2_max));     SaO2Ast = min(0.95, ((CO2_01    + CO2_12Ast)/2)./(CO2_max));
            ScO2NO = min(0.95, ((CO2_12NO + CO2_23NO)/2)./(CO2_max));     ScO2NPY = min(0.95, ((CO2_12NPY + CO2_23NPY)/2)./(CO2_max));     ScO2Pyr = min(0.95, ((CO2_12Pyr + CO2_23Pyr)/2)./(CO2_max));     ScO2Ast = min(0.95, ((CO2_12Ast + CO2_23Ast)/2)./(CO2_max));
            SvO2NO = min(0.95, ((CO2_23NO + CO2_34NO)/2)./(CO2_max));     SvO2NPY = min(0.95, ((CO2_23NPY + CO2_34NPY)/2)./(CO2_max));     SvO2Pyr = min(0.95, ((CO2_23Pyr + CO2_34Pyr)/2)./(CO2_max));     SvO2Ast = min(0.95, ((CO2_23Ast + CO2_34Ast)/2)./(CO2_max));

            % calculate the HbR
            HbRaNO = V1NO .* (1-SaO2NO);     HbRaNPY = V1NPY .* (1-SaO2NPY);        HbRaPyr = V1Pyr .* (1-SaO2Pyr);      HbRaAst = V1(1) .* (1-SaO2Ast);           
            HbRcNO = V2NO .* (1-ScO2NO);     HbRcNPY = V2NPY .* (1-ScO2NPY);        HbRcPyr = V2Pyr .* (1-ScO2Pyr);      HbRcAst = V2(1) .* (1-ScO2Ast);   
            HbRvNO = V3NO .* (1-SvO2NO);     HbRvNPY = V3NPY .* (1-SvO2NPY);        HbRvPyr = V3Pyr .* (1-SvO2Pyr);      HbRvAst = V3(1) .* (1-SvO2Ast);
        
            HbRNO  = HbRaNO  + HbRcNO  + HbRvNO ;
            HbRNPY = HbRaNPY + HbRcNPY + HbRvNPY;
            HbRPyr = HbRaPyr + HbRcPyr + HbRvPyr;
            HbRAst = HbRaAst + HbRcAst + HbRvAst;

            % due to the pysiological cap (95%) correct for lost HbR
            HbR  = simulationBOLD.(Experiment).y(:,idxHbR );
            missingHbRcontent = (HbR - HbR(1)) - (HbRNO - HbRNO(1) + HbRNPY - HbRNPY(1) + HbRPyr - HbRPyr(1) + HbRAst - HbRAst(1));
            
            % lost HbR is distributed based on the metabolic activity
            HbRNO  = HbRNO  + (abs(CMRO2NO )./CMRO2tot).*missingHbRcontent;
            HbRNPY = HbRNPY + (abs(CMRO2NPY)./CMRO2tot).*missingHbRcontent;
            HbRPyr = HbRPyr + (abs(CMRO2Pyr)./CMRO2tot).*missingHbRcontent;
            HbRAst = HbRAst + (abs(CMRO2Ast)./CMRO2tot).*missingHbRcontent;

            %% BOLD
            CBVNO  = V1NO  + V2NO  + V3NO;
            CBVNPY = V1NPY + V2NPY + V3NPY;
            CBVPyr = V1Pyr + V2Pyr + V3Pyr;
            
            BOLD_OISNO  = exp(-k_yBOLD*(HbRNO -HbRNO(1) ) - k_yBOLD2*(CBVNO  - CBVss));
            BOLD_OISNPY = exp(-k_yBOLD*(HbRNPY-HbRNPY(1)) - k_yBOLD2*(CBVNPY - CBVss));
            BOLD_OISPyr = exp(-k_yBOLD*(HbRPyr-HbRPyr(1)) - k_yBOLD2*(CBVPyr - CBVss));
            BOLD_OISAst = exp(-k_yBOLD*(HbRAst-HbRAst(1)) - k_yBOLD2*(CBVss  - CBVss));
           
            % calulate the area of the signal
            BOLD_OISNO_effect  =  trapz( abs( BOLD_OISNO(2:end ) -1));
            BOLD_OISNPY_effect =  trapz( abs( BOLD_OISNPY(2:end) -1));
            BOLD_OISPyr_effect =  trapz( abs( BOLD_OISPyr(2:end) -1));
            BOLD_OISAst_effect =  trapz( abs( BOLD_OISAst(2:end) -1));

            % calculate the percentual contribution of the total
            OISBOLD_NO_percent  = round( 100*( BOLD_OISNO_effect  / (BOLD_OISNO_effect + BOLD_OISNPY_effect + BOLD_OISPyr_effect + BOLD_OISAst_effect) ), 2);
            OISBOLD_NPY_percent = round( 100*( BOLD_OISNPY_effect / (BOLD_OISNO_effect + BOLD_OISNPY_effect + BOLD_OISPyr_effect + BOLD_OISAst_effect) ), 2);
            OISBOLD_Pyr_percent = round( 100*( BOLD_OISPyr_effect / (BOLD_OISNO_effect + BOLD_OISNPY_effect + BOLD_OISPyr_effect + BOLD_OISAst_effect) ), 2);
            OISBOLD_ASt_percent = round( 100*( BOLD_OISAst_effect / (BOLD_OISNO_effect + BOLD_OISNPY_effect + BOLD_OISPyr_effect + BOLD_OISAst_effect) ), 2);

            % assign value to cell array 
            BOLDTable{i,2}{idParam,1} = OISBOLD_NO_percent ;
            BOLDTable{i,3}{idParam,1} = OISBOLD_NPY_percent;
            BOLDTable{i,4}{idParam,1} = OISBOLD_Pyr_percent;
            BOLDTable{i,5}{idParam,1} = OISBOLD_ASt_percent;
            
            %% CBV effect
            CBVNO_effect  =  trapz( abs( CBVNO(2:end ) - CBVNO(2) ));
            CBVNPY_effect =  trapz( abs( CBVNPY(2:end) - CBVNPY(2)));
            CBVPyr_effect =  trapz( abs( CBVPyr(2:end) - CBVPyr(2)));
            CBVAst_effect =  0;
            
            % calculate the percentual contribution of the total
            CBV_NO_percent  = round( 100*( CBVNO_effect  / (CBVNO_effect + CBVNPY_effect + CBVPyr_effect + CBVAst_effect) ), 2);
            CBV_NPY_percent = round( 100*( CBVNPY_effect / (CBVNO_effect + CBVNPY_effect + CBVPyr_effect + CBVAst_effect) ), 2);
            CBV_Pyr_percent = round( 100*( CBVPyr_effect / (CBVNO_effect + CBVNPY_effect + CBVPyr_effect + CBVAst_effect) ), 2);
            CBV_ASt_percent = round( 100*( CBVAst_effect / (CBVNO_effect + CBVNPY_effect + CBVPyr_effect + CBVAst_effect) ), 2);
            
            % assign value to cell array 
            CBVTable{i,2}{idParam,1} = CBV_NO_percent ;
            CBVTable{i,3}{idParam,1} = CBV_NPY_percent;
            CBVTable{i,4}{idParam,1} = CBV_Pyr_percent;
            CBVTable{i,5}{idParam,1} = CBV_ASt_percent;

            %% CMRO2 effect
            CMRO2NO_effect  =  trapz( abs( CMRO2NO(2:end ) - CMRO2NO(2) ));
            CMRO2NPY_effect =  trapz( abs( CMRO2NPY(2:end) - CMRO2NPY(2)));
            CMRO2Pyr_effect =  trapz( abs( CMRO2Pyr(2:end) - CMRO2Pyr(2)));
            CMRO2Ast_effect =  trapz( abs( CMRO2Ast(2:end) - CMRO2Ast(2)));
            
            % calculate the percentual contribution of the total
            CMRO2_NO_percent  = round( 100*( CMRO2NO_effect  / (CMRO2NO_effect + CMRO2NPY_effect + CMRO2Pyr_effect + CMRO2Ast_effect) ), 2);
            CMRO2_NPY_percent = round( 100*( CMRO2NPY_effect / (CMRO2NO_effect + CMRO2NPY_effect + CMRO2Pyr_effect + CMRO2Ast_effect) ), 2);
            CMRO2_Pyr_percent = round( 100*( CMRO2Pyr_effect / (CMRO2NO_effect + CMRO2NPY_effect + CMRO2Pyr_effect + CMRO2Ast_effect) ), 2);
            CMRO2_ASt_percent = round( 100*( CMRO2Ast_effect / (CMRO2NO_effect + CMRO2NPY_effect + CMRO2Pyr_effect + CMRO2Ast_effect) ), 2);
            
            % assign value to cell array 
            CMRO2Table{i,2}{idParam,1} = CMRO2_NO_percent ;
            CMRO2Table{i,3}{idParam,1} = CMRO2_NPY_percent;
            CMRO2Table{i,4}{idParam,1} = CMRO2_Pyr_percent;
            CMRO2Table{i,5}{idParam,1} = CMRO2_ASt_percent;

            %% HbR effect
            HbRNO_effect  =  trapz( abs( HbRNO(2:end ) - HbRNO(2) ));
            HbRNPY_effect =  trapz( abs( HbRNPY(2:end) - HbRNPY(2)));
            HbRPyr_effect =  trapz( abs( HbRPyr(2:end) - HbRPyr(2)));
            HbRAst_effect =  trapz( abs( HbRAst(2:end) - HbRAst(2)));
            
            % calculate the percentual contribution of the total
            HbR_NO_percent  = round( 100*( HbRNO_effect  / (HbRNO_effect + HbRNPY_effect + HbRPyr_effect + HbRAst_effect) ), 2);
            HbR_NPY_percent = round( 100*( HbRNPY_effect / (HbRNO_effect + HbRNPY_effect + HbRPyr_effect + HbRAst_effect) ), 2);
            HbR_Pyr_percent = round( 100*( HbRPyr_effect / (HbRNO_effect + HbRNPY_effect + HbRPyr_effect + HbRAst_effect) ), 2);
            HbR_ASt_percent = round( 100*( HbRAst_effect / (HbRNO_effect + HbRNPY_effect + HbRPyr_effect + HbRAst_effect) ), 2);
            
            % assign value to cell array 
            HbRTable{i,2}{idParam,1} = HbR_NO_percent ;
            HbRTable{i,3}{idParam,1} = HbR_NPY_percent;
            HbRTable{i,4}{idParam,1} = HbR_Pyr_percent;
            HbRTable{i,5}{idParam,1} = HbR_ASt_percent;
            
            %% plot all experiments
            if and(shouldIPlot, isequal(idParam,1)) % plot the model behaviors for the best parameters
                % Volumes
                figure(1)
                sgtitle('V1')
                nexttile
                hold on
                plot(simulationBOLD.(Experiment).t, V1NO,    '-g', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, V1NPY,   '-r', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, V1Pyr,   '-b', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, sum([V1NO, V1NPY, V1Pyr],2) - 2*V1(1), 'color',  [.5 .5 .5], 'LineWidth', 4)
                plot(simulationBOLD.(Experiment).t, V1,      '-k', 'LineWidth', 1)
            
                xlabel('Time (s)');     ylabel('Vaso fractional');  title(Experiment, 'Interpreter', 'None');
                l = legend('V1NO', 'V1NPY', 'V1Pyr', 'Total', 'True');  l.EdgeColor = 'none';
                
                % CMRO2
                figure(2)
                sgtitle('CMRO_2')
                nexttile
                hold on
                plot(simulationBOLD.(Experiment).t, CMRO2(1) + CMRO2NO ,  '-g', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, CMRO2(1) + CMRO2NPY,  '-r', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, CMRO2(1) + CMRO2Pyr,  '-b', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, CMRO2(1) + CMRO2Ast,  '-m', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, CMRO2(1) + CMRO2NO + CMRO2NPY + CMRO2Pyr + CMRO2Ast, 'color', [.5 .5 .5], 'LineWidth', 4)
                plot(simulationBOLD.(Experiment).t, simulationBOLD.(Experiment).y(:,idxCMRO2), '-k', 'LineWidth', 1)
            
                xlabel('Time (s)');     ylabel('CMRO_2');   title(Experiment, 'Interpreter', 'None');
                l = legend('NO', 'NPY', 'Pyr', 'Astro', 'Total', 'True');   l.EdgeColor = 'none';
                
                % HbR
                figure(3)
                sgtitle('HbR')
                nexttile
                hold on
                plot(simulationBOLD.(Experiment).t, HbRNO  ,  '-g', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, HbRNPY ,  '-r', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, HbRPyr ,  '-b', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, HbRAst ,  '-m', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, sum([HbRNO HbRNPY HbRPyr HbRAst],2) - 3*simulationBOLD.(Experiment).y(1,idxHbR), 'Color', [.7 .7 .7], 'LineWidth', 4)
                plot(simulationBOLD.(Experiment).t, simulationBOLD.(Experiment).y(:,idxHbR)    ,  '-k', 'LineWidth', 1)
    
                xlabel('Time (s)');     ylabel('HbR');  title(Experiment, 'Interpreter', 'None')
                l = legend('NO', 'NPY', 'Pyr', 'Astro', 'Total', 'True');     l.EdgeColor = 'none';
                
                % BOLD
                figure(4)
                sgtitle('BOLD')
                nexttile
                hold on
                plot(simulationBOLD.(Experiment).t, BOLD_OISNO  , '-g', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, BOLD_OISNPY , '-r', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, BOLD_OISPyr , '-b', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, BOLD_OISAst , '-m', 'LineWidth', 2)
                plot(simulationBOLD.(Experiment).t, BOLD_OISNO .* BOLD_OISNPY .* BOLD_OISPyr .* BOLD_OISAst, 'Color', [.7 .7 .7], 'LineWidth', 4)
                plot(simulationBOLD.(Experiment).t, simulationBOLD.(Experiment).y(:,idxBOLD_OIS),  '-k', 'LineWidth', 1)
    
                xlabel('Time (s)');     ylabel('BOLD-OIS');     title(Experiment, 'Interpreter', 'None')
                l = legend('NO', 'NPY', 'Pyr', 'Astro', 'Total', 'True');     l.EdgeColor = 'none';
            end
        end
    end
end

save('Results/CellSpecificContributions/BOLDTable.mat' , 'BOLDTable' )
save('Results/CellSpecificContributions/CBVTable.mat'  , 'CBVTable'  )
save('Results/CellSpecificContributions/CMRO2Table.mat', 'CMRO2Table')
save('Results/CellSpecificContributions/HbRTable.mat'  , 'HbRTable'  )

%% PLot uncertainty
plotCBVcontributionFromFile()   % call plotfunction 
plotCMRO2contributionFromFile() % call plotfunction 
plotHbRcontributionFromFile()   % call plotfunction 
plotBOLDcontributionFromFile()  % call plotfunction 