clc
clear all

%PM
L0 = 18e-3;
S = 64e-6;
m = 6e-3;
epsilon = 0.95e+3;
L1 = 93e-3;
kappa = 0.25e-3;

%AK
% L0 = 39e-3;
% S = 50e-6;
% m = 8e-3;
% epsilon = 6.65e+3;
% L1 = 415e-3;
% kappa = 1.75e-3;

gamma = 1.3;
patm = 1e+5;

deltaU = epsilon * (1 - (L0/L1)^(gamma-1))
K = deltaU
KPD = K / epsilon * 100
v = sqrt(2*K/m)
p0 = (gamma-1)*epsilon/S/L0 / patm
p1 = p0*(L0/L1)^gamma
rho = 3/2 * kappa * 44/101 / (L0*S)
