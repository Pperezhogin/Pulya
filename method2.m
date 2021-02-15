clc
clear all

N = 64;
h = 1 / (N-1);
y = 0:h:1;

%PM
% L0 = 18e-3;
% S = 64e-6;
% m = 6e-3;
% epsilon = 0.95e+3;
% L1 = 93e-3;
% kappa = 0.25e-3;

% AK47
L0 = 39e-3;
S = 50e-6;
m = 8e-3;
epsilon = 6.65e+3;
L1 = 415e-3;
kappa = 1.75e-3;

rho0 = 3/2 * kappa * 44/101 / (L0*S);
gamma = 1.3;
patm = 1e+5;

rho0 = rho0;

K = epsilon * (1 - (L0/L1)^(gamma-1));
v_estimate = sqrt(2*K/m);

rho = rho0 * ones(size(y));
rhou = zeros(size(y));
rhoE = epsilon / S / L0 * ones(size(y));

t = 0;
L = L0;
v = 0;
CFL = 0.1;

vn(1) = 0;
tn(1) = 0;
n = 1;
while (L<L1)
    %[rho,rhou,rhoE,L,v,t,dt,u,p,e,c_sound,mass] = collocated_neuman_bc(rho, rhou, rhoE, m, S, L, v, gamma, CFL, t);
    [rho,rhou,rhoE,L,v,t,dt,u,p,e,c_sound,mass] = collocated_conservative_bc(rho, rhou, rhoE, m, S, L, v, gamma, CFL, t);
    %[rho,rhou,rhoE,L,v,t,dt,u,p,e,c_sound,mass] = staggered_conservative_bc(rho, rhou, rhoE, m, S, L, v, gamma, CFL, t);
    n = n + 1;
    vn(n) = v;
    tn(n) = t;
    mass
end

%plot(u)

plot(tn,vn)
xlabel('time ,sec')
ylabel('velocity, m/s')
title('The AK-47')
grid on
hold on
plot(tn,v_estimate*ones(size(tn)),'--','LineWidth',2)
