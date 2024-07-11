function [model] = Sten2023_v31_steadyState_syms()

% set the parametrisation of the problem options are 'log', 'log10' and

model.param = 'log10';

%%
% STATES
% create state syms 
syms N_NO N_NPY N_Pyr N_Astro N_SOM GABA Glut Ca_NO Ca_NPY Ca_Pyr Ca_Astro Ca_SOM AA PGE2 PGE2vsm SOM_NO SOM_NPY NO NOvsm NPY NPYvsm V1 V2 V3 f1 f2 f3 nO2_1 nO2_2 nO2_3 nO2_t
% create state vector
model.sym.x = [N_NO N_NPY N_Pyr N_Astro N_SOM GABA Glut Ca_NO Ca_NPY Ca_Pyr Ca_Astro Ca_SOM AA PGE2 PGE2vsm SOM_NO SOM_NPY NO NOvsm NPY NPYvsm V1 V2 V3 f1 f2 f3 nO2_1 nO2_2 nO2_3 nO2_t];
%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms 

syms k_u1 k_u2 k_u3 k_u4 N_NOMax N_NPYMax N_PyrMax N_SOMMax                                                                 % #1-8   : Stimulation and saturation parameters
syms kPF1 kPF2 kPF3 kPF4 kPF5 kIN kIN2 kIN3 kIN4 kIN5 sinkN_NO sinkN_NPY sinkN_Pyr sinkN_Astro sinkN_SOM                     % #9-23  : Electrical signalling parameters 
syms kGABA_NO kGABA_NPY kGABA_Astro kGABA_SOM sink_GABA kGlut_Astro kGlut_Pyr sink_Glut                                     % #24-31 : GABA and Glutatmate signalling  
syms k_Ca k_CaNO k_CaNPY k_CaPyr k_CaAstro kPFCa kPF2Ca k_CaSOM sinkCa_NO sinkCa_NPY sinkCa_Pyr sinkCa_Astro sinkCa_SOM     % #32-44 : Ca dynamics parameters
syms kPL kCOX Km_AA kPGE2 Km_PGE2 sinkPGE2 kNOS kNO Km_NO kNOS2 kNO2 Km_NO2 sinkNO kNPY Vmax_NPY Km_NPY kNPY2 Vmax_NPY2 Km_NPY2 sinkNPY   % #45-64 : Intracellular signalling parameters
syms ky1 ky2 ky3 K1 K2 K3 vis1 vis2 vis3 k_leak k_capperm                                                                   % #65-75 : Vascular parameters
syms kCMRO2_NO kCMRO2_NPY kCMRO2_Pyr kCMRO2_Astro kCMRO2_SOM k_yBOLD k_yBOLD2                                             % #76-82 : CMRO2 and BOLD-scaling parameters

% create parameter vector 
model.sym.p = [k_u1 k_u2 k_u3 k_u4 N_NOMax N_NPYMax N_PyrMax N_SOMMax kPF1 kPF2 kPF3 kPF4 kPF5 kIN kIN2 kIN3 kIN4 kIN5 sinkN_NO sinkN_NPY sinkN_Pyr sinkN_Astro sinkN_SOM kGABA_NO kGABA_NPY kGABA_Astro kGABA_SOM sink_GABA kGlut_Astro kGlut_Pyr sink_Glut k_Ca k_CaNO k_CaNPY k_CaPyr k_CaAstro kPFCa kPF2Ca k_CaSOM sinkCa_NO sinkCa_NPY sinkCa_Pyr sinkCa_Astro sinkCa_SOM kPL kCOX Km_AA kPGE2 Km_PGE2 sinkPGE2 kNOS kNO Km_NO kNOS2 kNO2 Km_NO2 sinkNO kNPY Vmax_NPY Km_NPY kNPY2 Vmax_NPY2 Km_NPY2 sinkNPY ky1 ky2 ky3 K1 K2 K3 vis1 vis2 vis3 k_leak k_capperm kCMRO2_NO kCMRO2_NPY kCMRO2_Pyr kCMRO2_Astro kCMRO2_SOM k_yBOLD k_yBOLD2 ];
%%  
% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create constants syms
syms u1 u2 u3 u4 NOvsm0 PGE2vsm0 NPYvsm0 kCa g_1SS g_2SS g_3SS g_s CMRO2_0 CO2_lSS pO2_femart HbO_0 HbR_0 SaO2_0 ScO2_0 SvO2_0 TE B0  % kCa is unused 
% create constants vector 
model.sym.k = [u1 u2 u3 u4 NOvsm0 PGE2vsm0 NPYvsm0 kCa g_1SS g_2SS g_3SS g_s CMRO2_0 CO2_lSS pO2_femart HbO_0 HbR_0 SaO2_0 ScO2_0 SvO2_0 TE B0];

%%
% SYSTEM EQUATIONS
% create symbolic variable for time
syms t 
model.sym.xdot = sym(zeros(size(model.sym.x)));

model.sym.xdot(1) =  + kPF1*Glut - kIN*GABA  - sinkN_NO*N_NO;        % N_NO
model.sym.xdot(2) =  + kPF2*Glut - kIN2*GABA - sinkN_NPY*N_NPY;      % N_NPY
model.sym.xdot(3) =  + kPF3*Glut - kIN3*GABA - sinkN_Pyr*N_Pyr;      % N_Pyr
model.sym.xdot(4) =  + kPF4*Glut - kIN4*GABA - sinkN_Astro*N_Astro;  % N_Astro
model.sym.xdot(5) =  + kPF5*Glut - kIN5*GABA - sinkN_SOM*N_SOM;      % N_SOM

model.sym.xdot(6) = + kGABA_NO*Ca_NO + kGABA_NPY*Ca_NPY + kGABA_Astro*Ca_Astro + kGABA_SOM*Ca_SOM - sink_GABA*GABA;
model.sym.xdot(7) = + kGlut_Astro*Ca_Astro + kGlut_Pyr*Ca_Pyr - sink_Glut*Glut;

%% Calcium dynamics
model.sym.xdot(8) =  log(1 + exp(k_Ca + k_CaNO*N_NO))                                  - sinkCa_NO*Ca_NO;       % Ca_NO
model.sym.xdot(9) =  log(1 + exp(k_Ca + k_CaNPY*N_NPY))                                - sinkCa_NPY*Ca_NPY;     % Ca_NPY
model.sym.xdot(10) = log(1 + exp(k_Ca + k_CaPyr*N_Pyr))                                - sinkCa_Pyr*Ca_Pyr;     % Ca_Pyr
model.sym.xdot(11) = log(1 + exp(k_Ca + k_CaAstro*N_Astro)) + kPFCa*Glut + kPF2Ca*GABA - sinkCa_Astro*Ca_Astro; % Ca_Astro
model.sym.xdot(12) = log(1 + exp(k_Ca + k_CaSOM*N_SOM))                                - sinkCa_SOM*Ca_SOM;     % Ca_SOM

%% Pyramidal intracellular signaling AA->PGE2
model.sym.xdot(13) = kPL*Ca_Pyr - kCOX*AA/(Km_AA+AA);    % AA
model.sym.xdot(14) = kCOX*AA/(Km_AA+AA) - kPGE2*(PGE2/(Km_PGE2 + PGE2));    % PGE2
model.sym.xdot(15) = kPGE2*(PGE2/(Km_PGE2 + PGE2)) - sinkPGE2*PGE2vsm;      % PGE2vsm

%% GABAergic SOM signaling
model.sym.xdot(16) = kNOS2*Ca_SOM - kNO2     *(SOM_NO^1/(Km_NO2^1   + SOM_NO^1));% SOM_NO
model.sym.xdot(17) = kNPY2*Ca_SOM - Vmax_NPY2* SOM_NPY /(Km_NPY2 + SOM_NPY);     % SOM_NPY

%% GABAergic NO signaling
model.sym.xdot(18) = kNOS*Ca_NO - kNO*(NO^1/(Km_NO^1 + NO^1));     % NO
model.sym.xdot(19) = kNO*(NO^1/(Km_NO^1 + NO^1)) + kNO2*(SOM_NO^1/(Km_NO2^1 + SOM_NO^1)) - sinkNO*NOvsm;   % NOvsm

%% GABAergic NPY signaling
model.sym.xdot(20) = kNPY*Ca_NPY - Vmax_NPY*NPY/(Km_NPY+NPY);     % NPY
model.sym.xdot(21) = Vmax_NPY*NPY/(Km_NPY+NPY) + Vmax_NPY2* SOM_NPY/(Km_NPY2 + SOM_NPY) - sinkNPY*NPYvsm;  % NPYvsm

%% Blood Vessel circuit
% V1ss = 0.29;     V2ss = 0.44;     V3ss = 0.27; % steady-state volume
% R1ss = 0.74;     R2ss = 0.08;     R3ss = 0.18; % steay-state resistance
fSS = 1; 

model.sym.xdot(22) = 0; % V1 - Arterial CBV
model.sym.xdot(23) = 0; % V2 - Capillary CBV
model.sym.xdot(24) = 0; % V3 - Venous CBV

f0 = 1;
model.sym.xdot(25) = 0; % f1 - Arterial CBF
model.sym.xdot(26) = 0; % f2 - Capillary CBF
model.sym.xdot(27) = 0; % f3 - Venous CBF


%% Barret Hb and blood pressure 
CO2_max = 9.26;
V_t = 34.8;
p50 = 36;           h = 2.6;
sigma_O2 = 1.46*10^-3;   

% dynamic O2 leakage function
CO2_l = am_max(CO2_lSS*(1 - k_leak*(f0 - fSS)), 0); %CO2 leak is not allowed to be negative

% Concentrations of O2 
CO2_0 = CO2_max/(10^(log10(pO2_femart/p50)/(-1/h))+1);
CO2_01 = CO2_0 - CO2_l;
CO2_12 = nO2_1/V1;
CO2_23 = nO2_2/V2;
CO2_34 = nO2_3/V3;
CO2_t = nO2_t/V_t;

% partial oxygen preassure 
pO2_01 = p50*((CO2_max/CO2_01) -1)^(-1/h);
pO2_12 = p50*((CO2_max/CO2_12) -1)^(-1/h);
pO2_23 = p50*((CO2_max/CO2_23) -1)^(-1/h);  
pO2_34 = p50*((CO2_max/CO2_34) -1)^(-1/h);

%oxygen preassure
pO2_1 = (pO2_01 + pO2_12)/2;
pO2_2 = (pO2_12 + pO2_23)/2;
pO2_3 = (pO2_23 + pO2_34)/2;
pO2_t = CO2_t/sigma_O2;


g_2 = g_2SS*(1 + k_capperm*(f0-fSS)); % dynamic capilary permability function  

% differance in patial oxygen preassure
jO2_1 = g_1SS*(pO2_1 - pO2_t);
jO2_2 = g_2*(pO2_2 - pO2_t);
jO2_3 = g_3SS*(pO2_3 - pO2_t);
jO2_s = g_s*(pO2_1 - pO2_3);

%Metabolic rate of oxyen
CMRO2 = CMRO2_0;

model.sym.xdot(28) = f0*CO2_01 - f1*CO2_12 - jO2_1 - jO2_s; % nO2_1 - amount of Arterial oxygen   
model.sym.xdot(29) = f1*CO2_12 - f2*CO2_23 - jO2_2;         % nO2_2 - amount of Capillary oxygen
model.sym.xdot(30) = f2*CO2_23 - f3*CO2_34 - jO2_3 + jO2_s; % nO2_3 - amount of venous oxygen
model.sym.xdot(31) = jO2_1 + jO2_2 + jO2_3 - CMRO2;         % nO2_t - total amount of O2

%oxygen saturatuion
SaO2 = ((CO2_01 + CO2_12)/2)/(CO2_max);
ScO2 = ((CO2_12 + CO2_23)/2)/(CO2_max); 
SvO2 = ((CO2_23 + CO2_34)/2)/(CO2_max);

% Hemoglobin per vessel compartment Oxygenated and deoxygenated 
HbOa = V1*SaO2;     HbRa = V1*(1-SaO2);          
HbOc = V2*ScO2;     HbRc = V2*(1-ScO2);           
HbOv = V3*SvO2;     HbRv = V3*(1-SvO2);        

HbO = HbOa + HbOc + HbOv;
HbR = HbRa + HbRc + HbRv;

%%
% INITIAL CONDITIONS
model.sym.x0(1:21) = 0;
model.sym.x0(22:24) = [0.29; 0.44; 0.27];
model.sym.x0(25:27) = 1;
model.sym.x0(28:31) = 1;

model.sym.dx0 = sym(zeros(size(model.sym.x)));

%%
% OBSERVALES
syms OxyHb DeOxyHb artO2Sat capO2Sat venO2sat
model.sym.y = [OxyHb DeOxyHb artO2Sat capO2Sat venO2sat]; 
 
model.sym.y(1) = HbO;
model.sym.y(2) = HbR;
model.sym.y(3) = SaO2;
model.sym.y(4) = ScO2;
model.sym.y(5) = SvO2;
end