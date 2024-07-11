function [model] = Sten2023_v21_BOLDcontributions_syms()

% set the parametrisation of the problem options are 'log', 'log10' and

model.param = 'log10';

%%
% STATES
% create state syms 
syms N_NO N_NPY N_Pyr N_Astro GABA Glut Ca_NO Ca_NPY Ca_Pyr Ca_Astro AA PGE2 PGE2vsm NO NOvsm NPY NPYvsm V1 V2 V3 f1 f2 f3 nO2_1 nO2_2 nO2_3 nO2_t
% create state vector
model.sym.x = [N_NO N_NPY N_Pyr N_Astro GABA Glut Ca_NO Ca_NPY Ca_Pyr Ca_Astro AA PGE2 PGE2vsm NO NOvsm NPY NPYvsm V1 V2 V3 f1 f2 f3 nO2_1 nO2_2 nO2_3 nO2_t];

%%
% INITIAL CONDITIONS
model.sym.x0(1:17) = 0;
model.sym.x0(18:20) = [0.29; 0.44; 0.27];
model.sym.x0(21:23) = 1;
model.sym.x0(24:27) = 1;

model.sym.dx0 = sym(zeros(size(model.sym.x)));
%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms 
syms k_u1 k_u2 k_u3 N_NOMax N_NPYMax N_PyrMax                                                           % #1-6   : Stimulation and saturation parameters
syms kPF1 kPF2 kPF3 kPF4 kIN kIN2 kIN3 kIN4 sinkN_NO sinkN_NPY sinkN_Pyr sinkN_Astro                    % #7-18  : Electrical signalling parameters 
syms kGABA_NO kGABA_NPY kGABA_Astro sink_GABA kGlut_Astro kGlut_Pyr sink_Glut                           % #19-25 : GABA and Glutatmate signalling  
syms k_Ca k_CaNO k_CaNPY k_CaPyr k_CaAstro kPFCa kPF2Ca sinkCa_NO sinkCa_NPY sinkCa_Pyr sinkCa_Astro    % #26-36 : Ca dynamics parameters
syms kPL kCOX Km_AA kPGE2 Km_PGE2 sinkPGE2 kNOS kNO Km_NO sinkNO kNPY Vmax_NPY Km_NPY sinkNPY           % #37-50 : Intracellular signalling parameters
syms ky1 ky2 ky3 K1 K2 K3 vis1 vis2 vis3 k_leak k_capperm                                               % #51-61 : Vascular parameters
syms kCMRO2_NO kCMRO2_NPY kCMRO2_Pyr kCMRO2_Astro k_yBOLD k_yBOLD2                                      % #62-67 : CMRO2 and BOLD-scaling parameters

% create parameter vector 
model.sym.p = [k_u1 k_u2 k_u3 N_NOMax N_NPYMax N_PyrMax kPF1 kPF2 kPF3 kPF4 kIN kIN2 kIN3 kIN4 sinkN_NO sinkN_NPY sinkN_Pyr sinkN_Astro kGABA_NO kGABA_NPY kGABA_Astro sink_GABA kGlut_Astro kGlut_Pyr sink_Glut k_Ca k_CaNO k_CaNPY k_CaPyr k_CaAstro kPFCa kPF2Ca sinkCa_NO sinkCa_NPY sinkCa_Pyr sinkCa_Astro kPL kCOX Km_AA kPGE2 Km_PGE2 sinkPGE2 kNOS kNO Km_NO sinkNO kNPY Vmax_NPY Km_NPY sinkNPY ky1 ky2 ky3 K1 K2 K3 vis1 vis2 vis3 k_leak k_capperm kCMRO2_NO kCMRO2_NPY kCMRO2_Pyr kCMRO2_Astro k_yBOLD k_yBOLD2];
%%  
% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
syms u1 u2 u3 NOvsm0 PGE2vsm0 NPYvsm0 kCa g_1SS g_2SS g_3SS g_s CMRO2_0 CO2_lSS pO2_femart HbO_0 HbR_0 SaO2_0 ScO2_0 SvO2_0 TE B0  % kCa is unused 
% create parameter vector 
model.sym.k = [u1 u2 u3 NOvsm0 PGE2vsm0 NPYvsm0 kCa g_1SS g_2SS g_3SS g_s CMRO2_0 CO2_lSS pO2_femart HbO_0 HbR_0 SaO2_0 ScO2_0 SvO2_0 TE B0];

%%
% SYSTEM EQUATIONS
% create symbolic variable for time
syms t 
model.sym.xdot = sym(zeros(size(model.sym.x)));

model.sym.xdot(1) =  + k_u1*u1*am_max(N_NOMax-N_NO,0)   + kPF1*Glut - kIN*GABA  - sinkN_NO*N_NO;        % N_NO
model.sym.xdot(2) =  + k_u2*u2*am_max(N_NPYMax-N_NPY,0) + kPF2*Glut - kIN2*GABA - sinkN_NPY*N_NPY;      % N_NPY
model.sym.xdot(3) =  + k_u3*u3*am_max(N_PyrMax-N_Pyr,0) + kPF3*Glut - kIN3*GABA - sinkN_Pyr*N_Pyr;      % N_Pyr
model.sym.xdot(4) =                                     + kPF4*Glut - kIN4*GABA - sinkN_Astro*N_Astro;  % N_Astro

model.sym.xdot(5) = + kGABA_NO*Ca_NO + kGABA_NPY*Ca_NPY + kGABA_Astro*Ca_Astro - sink_GABA*GABA;
model.sym.xdot(6) = + kGlut_Astro*Ca_Astro + kGlut_Pyr*Ca_Pyr - sink_Glut*Glut;

%% Calcium dynamics
model.sym.xdot(7) =  log(1 + exp(k_Ca + k_CaNO*N_NO))                                  - sinkCa_NO*Ca_NO;       % Ca_NO
model.sym.xdot(8) =  log(1 + exp(k_Ca + k_CaNPY*N_NPY))                                - sinkCa_NPY*Ca_NPY;     % Ca_NPY
model.sym.xdot(9) =  log(1 + exp(k_Ca + k_CaPyr*N_Pyr))                                - sinkCa_Pyr*Ca_Pyr;     % Ca_Pyr
model.sym.xdot(10) = log(1 + exp(k_Ca + k_CaAstro*N_Astro)) + kPFCa*Glut + kPF2Ca*GABA - sinkCa_Astro*Ca_Astro; % Ca_Astro

%% Pyramidal intracellular signaling AA->PGE2
model.sym.xdot(11) = kPL*Ca_Pyr - kCOX*AA/(Km_AA+AA);    % AA
model.sym.xdot(12) = kCOX*AA/(Km_AA+AA) - kPGE2*(PGE2/(Km_PGE2 + PGE2));    % PGE2
model.sym.xdot(13) = kPGE2*(PGE2/(Km_PGE2 + PGE2)) - sinkPGE2*PGE2vsm;      % PGE2vsm

%% GABAergic NO signaling
model.sym.xdot(14) = kNOS*Ca_NO - kNO*(NO/(Km_NO + NO));     % NO
model.sym.xdot(15) = kNO*(NO/(Km_NO + NO)) - sinkNO*NOvsm;   % NOvsm

%% GABAergic NPY signaling
model.sym.xdot(16) = kNPY*Ca_NPY - Vmax_NPY*NPY/(Km_NPY+NPY);     % NPY
model.sym.xdot(17) = Vmax_NPY*NPY/(Km_NPY+NPY) - sinkNPY*NPYvsm;  % NPYvsm

stim_circNO  =  ky1*(NOvsm-NOvsm0);
stim_circNPY = -ky3*(NPYvsm-NPYvsm0);
stim_circPyr =  ky2*(PGE2vsm-PGE2vsm0);
stim_circ = stim_circNO + stim_circPyr + stim_circNPY;

%% Blood Vessel circuit
% 1 - Arterioles, 2 - Capillaries, 3 - Venules  
syms C1 C2 C3 f0 R1 R2 R3

V1ss = 0.29;     V2ss = 0.44;     V3ss = 0.27; % steady-state volume
R1ss = 0.74;     R2ss = 0.08;     R3ss = 0.18; % steay-state resistance
fSS = 1;

p1=1;                           p2=1-R1ss;                      p3=1-(R1ss+R2ss);
C1ss = V1ss/(p1 - 0.5*R1ss);    C2ss = V2ss/(p2 - 0.5*R2ss);    C3ss = V3ss/(p3 - 0.5*R3ss);% Vessel compliance 
L1=(R1ss*(V1ss^2))^(1/3);       L2=(R2ss*(V2ss^2))^(1/3);       L3=(R3ss*(V3ss^2))^(1/3);
R1=(L1^3)/(V1^2);               R2=(L2^3)/(V2^2);               R3=(L3^3)/(V3^2);           % Vessel resistance 

f0 = (2/R1)-((R1+R2)*f1 + (R2+R3)*f2 +R3*f3)/R1; % basal blood flow

% preassures
P1 = (1/2)*(R1 + R2)*f1 + (1/2)*(R2 + R3)*f2 + (1/2)*R3*f3; % P1 - preassure at C1
P2 = (1/2)*(R2 + R3)*f2 + (1/2)*R3*f3;                      % P2 - preassure at C2
P3 = (1/2)*R3*f3;                                           % P3 - preassure at C3

delta_p1 = (f0 + f1)*R1/2; % preassure drop over Arterial compartment 
delta_p2 = (f1 + f2)*R2/2; % preassure drop over capillary compartment 
delta_p3 = (f2 + f3)*R3/2; % preassure drop over venular compartment 

% Vessel compliances
C1 = V1/((1/2)*(R1+R2)*f1 + (1/2)*(R2+R3)*f2 + (R3/2)*f3);
C2 = V2/(                   (1/2)*(R2+R3)*f2 + (R3/2)*f3);
C3 = V3/(                                      (R3/2)*f3);

model.sym.xdot(18) = (1/vis1)*(((K1 -(V1/V1ss))/(K1-1)) + stim_circ - 2*(V1/(C1ss*((R1+R2)*f1+(R2+R3)*f2+R3*f3))));                 % V1 - Arterial CBV
model.sym.xdot(19) = (1/vis2)*(((K2 -(V2/V2ss))/(K2-1))             - 2*(V2/(C2ss*((R2+R3)*f2+R3*f3))));                            % V2 - Capillary CBV
model.sym.xdot(20) = (1/vis3)*(((K3 -(V3/V3ss))/(K3-1))             - 2*(V3/(C3ss*(R3*f3))));                                       % V3 - Venous CBV

model.sym.xdot(21) = f0 - model.sym.xdot(18) - f1;   % f1 - Arterial CBF
model.sym.xdot(22) = f1 - model.sym.xdot(19) - f2;   % f2 - Capillary CBF
model.sym.xdot(23) = f2 - model.sym.xdot(20) - f3;   % f3 - Venous CBF

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
jO2_2 = g_2*  (pO2_2 - pO2_t);
jO2_3 = g_3SS*(pO2_3 - pO2_t);
jO2_s = g_s*  (pO2_1 - pO2_3);

%Metabolic rate of oxyen
CMRO2 = CMRO2_0*(1 + kCMRO2_NO*(Ca_NO-model.sym.x0(7)) + kCMRO2_NPY*(Ca_NPY-model.sym.x0(8)) + kCMRO2_Pyr*(Ca_Pyr-model.sym.x0(9)) + kCMRO2_Astro*(Ca_Astro-model.sym.x0(10)));

model.sym.xdot(24) = f0*CO2_01 - f1*CO2_12 - jO2_1 - jO2_s; % nO2_1 - amount of Arterial oxygen   
model.sym.xdot(25) = f1*CO2_12 - f2*CO2_23 - jO2_2;         % nO2_2 - amount of Capillary oxygen
model.sym.xdot(26) = f2*CO2_23 - f3*CO2_34 - jO2_3 + jO2_s; % nO2_3 - amount of venous oxygen
model.sym.xdot(27) = jO2_1 + jO2_2 + jO2_3 - CMRO2;         % nO2_t - amount of tissue oxygen

%oxygen saturaturation
SaO2 = ((CO2_01 + CO2_12)/2)/(CO2_max);
ScO2 = ((CO2_12 + CO2_23)/2)/(CO2_max);
SvO2 = ((CO2_23 + CO2_34)/2)/(CO2_max);

% Hemoglobin per vessel compartment Oxygenated and deoxygenated 
HbOa = V1*SaO2;     HbRa = V1*(1-SaO2);          
HbOc = V2*ScO2;     HbRc = V2*(1-ScO2);           
HbOv = V3*SvO2;     HbRv = V3*(1-SvO2);        

HbO = HbOa + HbOc + HbOv;
HbR = HbRa + HbRc + HbRv;

%% BOLD Buxton  
VI = 0.05;
Va0 = V1ss*VI;      Vc0 = V2ss*VI;      Vv0 = V3ss*VI;
Va = V1*VI;         Vc = V2*VI;         Vv = V3*VI;
Ve = 1 - (Va+Vc+Vv);

% Constants 
Hct = 0.44;             Hct_c = 0.33;                   deltaChi = 2.64*10^-7;        
gamma = 2.68*10^8;      Cav = 302.06*Hct + 41.83;       Cc = 302.06*Hct_c + 41.83; 
Aav = 14.87*Hct+14.686; Ac = 14.87*Hct_c+14.686; 
R2e0 = 25.1;            R2a0 = Aav+Cav*((1-SaO2_0).^2); R2c0 = Ac+Cc*((1-ScO2_0).^2);   R2v0 = Aav+Cav*((1-SvO2_0).^2);

lambda = 1.15;  Se0 = exp(-TE*R2e0);    Sa0 = exp(-TE*R2a0);    Sc0 = exp(-TE*R2c0);    Sv0 = exp(-TE*R2v0);

epsA = lambda*(Sa0/Se0);    epsC = lambda*(Sc0/Se0); epsV = lambda*(Sv0/Se0);
%epsA=1.3;           epsC=1.02;          epsV=0.5;
preAV = 4*pi*Hct*deltaChi*gamma*B0/3;     preC = 0.04*((deltaChi*Hct_c*gamma*B0).^2);

% saturation for equal tissue–blood susceptibility
Yoff = 0.95;

%Vascular
deltaR2a = Cav*((1-SaO2).^2-((1-SaO2_0).^2));
deltaR2c = Cc*((1-ScO2).^2-((1-ScO2_0).^2));
deltaR2v = Cav*((1-SvO2).^2-((1-SvO2_0).^2));

% Extravascular R2
deltaR2ea = Va*(abs(Yoff-SaO2))-Va0*(abs(Yoff-SaO2_0));
deltaR2ev = Vv*(abs(Yoff-SvO2))-Vv0*(abs(Yoff-SvO2_0));
deltaR2ec = Vc*((abs(Yoff-ScO2)).^2)-Vc0*((abs(Yoff-ScO2_0)).^2);
deltaR2e  = preAV*(deltaR2ea+deltaR2ev) + preC*deltaR2ec;

%Signal from each compartment
Sa = epsA*Va*exp(-TE*deltaR2a);
Sc = epsC*Vc*exp(-TE*deltaR2c);
Sv = epsV*Vv*exp(-TE*deltaR2v);
Se = Ve*exp(-TE*deltaR2e);

H = ((1 - VI) + epsA*Va0 + epsC*Vc0 + epsV*Vv0);

%% M matrix
    matris = eye(size(model.sym.x,2),size(model.sym.x,2));
    matris(21,21) = 0;  matris(22,22) = 0; matris(23,23) = 0; 
    model.sym.M= matris;

%% output simplifications

CBV = V1 + V2 + V3;
CBVss = V1ss + V2ss + V3ss;
CBF = (V1*(f0 + f1) + V2*(f1 + f2) + V3*(f2 + f3))/(2*(V1 + V2 + V3));

CBFprecent = 100*(CBF-1);
CBVpercent = 100*(CBV-1);

BOLDfractional = (1/H)*(Se+Sa+Sc+Sv); % fractional BOLD change
BOLDpercent = 100*(BOLDfractional-1); % Precental BOLD change

BOLD_OIS = exp(-k_yBOLD*(HbR-HbR_0) - k_yBOLD2*(CBV - CBVss)); % proxy BOLD signal from deoxyHemoglobin
BOLD_OISpercent = 100*(BOLD_OIS -1);


HbTpercent = (((V1-V1ss) + (V2-V2ss) + (V3-V3ss))/(V1ss + V2ss + V3ss))*100; 
HbOpercent = ((HbO - HbO_0)/HbO_0)*100;                   
HbRpercent = ((HbR - HbR_0)/HbR_0)*100;                   

HbTrelative = 60*((HbO - HbO_0)/HbO_0) + 40*((HbR - HbR_0)/HbR_0); 
HbOrelative = 60*((HbO - HbO_0)/HbO_0);                   
HbRrelative = 40*((HbR - HbR_0)/HbR_0); 

%%
% OBSERVALES
syms           BOLD_OIS_y BOLD_OISpercent_y CBF_y CBV_y HbT_y HbO_y HbR_y HbTrelative_y HbOrelative_y HbRrelative_y BOLD_fMRI CMRO2_y HbR_y jO2_1_y jO2_2_y jO2_3_y CO2_01_y
model.sym.y = [BOLD_OIS_y BOLD_OISpercent_y CBF_y CBV_y HbT_y HbO_y HbR_y HbTrelative_y HbOrelative_y HbRrelative_y BOLD_fMRI CMRO2_y HbR_y jO2_1_y jO2_2_y jO2_3_y CO2_01_y]; 

model.sym.y(1 ) = BOLD_OIS; 
model.sym.y(2 ) = BOLD_OISpercent; 
model.sym.y(3 ) = CBFprecent; 
model.sym.y(4 ) = CBVpercent; % precental CBV
model.sym.y(5 ) = HbTpercent; % precental Hemoglobin Total 
model.sym.y(6 ) = HbOpercent; % precental OxyHemoglobin  
model.sym.y(7 ) = HbRpercent; % precental DeOxyHemoglobin 
model.sym.y(8 ) = HbTrelative;
model.sym.y(9 ) = HbOrelative;
model.sym.y(10) = HbRrelative;
model.sym.y(11) = BOLDpercent;   % Precental BOLD fMRI - signal

model.sym.y(12) = CMRO2; 
model.sym.y(13) = HbR; 
model.sym.y(14) = jO2_1; 
model.sym.y(15) = jO2_2; 
model.sym.y(16) = jO2_3; 
model.sym.y(17) = CO2_01; 

end