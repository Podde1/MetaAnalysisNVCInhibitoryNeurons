%% function to solve the intial values needed for the simulation of the circuit
function [g_1, g_2, g_3, g_s, CMRO2_0, CO2_l] = BarretPO2constants(pO2_01,pO2_12,pO2_23,pO2_34,PO2_t,pO2_femart)

% constants
CO2_max = 9.26;      
p50 = 36;           h = 2.6;

%% pressure unique values
CO2_0= CO2_max/(10^(log10(pO2_femart/p50)/(-1/h))+1);
CO2_01 = CO2_max/(10^(log10(pO2_01/p50)/(-1/h))+1);
CO2_12 = CO2_max/(10^(log10(pO2_12/p50)/(-1/h))+1);
CO2_23 = CO2_max/(10^(log10(pO2_23/p50)/(-1/h))+1);
CO2_34 = CO2_max/(10^(log10(pO2_34/p50)/(-1/h))+1);
CO2_l = CO2_0 - CO2_01;

PO2_1 = (pO2_01 + pO2_12)/2;
PO2_2 = (pO2_12 + pO2_23)/2;
PO2_3 = (pO2_23 + pO2_34)/2;

g_2 = (CO2_12 - CO2_23)/(PO2_2 - PO2_t);

% intersection of ekv 1 and ekv3 gives value for g_s
syms x y
eqn1 = (CO2_01 - CO2_12 - x*(PO2_1 - PO2_t))/(PO2_1 - PO2_3) == y;
eqn2 = (CO2_23 - CO2_34 - x*(PO2_3 - PO2_t))/(-1*(PO2_1 - PO2_3)) == y;

[A,B] = equationsToMatrix([eqn1,eqn2],[x,y]);

X = linsolve(A,B);
g_s = eval(X(1));

%now solve rest of g constants and CMRO2_0
g_1 = (CO2_01 - CO2_12 - g_s*(PO2_1 - PO2_3))/(PO2_1 - PO2_t);
g_3 = (CO2_23 - CO2_34 + g_s*(PO2_1 - PO2_3))/(PO2_3 - PO2_t);
CMRO2_0 = g_1*(PO2_1 - PO2_t) + g_2*(PO2_2 - PO2_t) + g_3*(PO2_3 - PO2_t);

end