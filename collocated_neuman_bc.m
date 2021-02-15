function [rho1,rhou1,rhoE1,L1,v1,t1,dt,u,p,e,c_sound,mass] = collocated_neuman_bc(rho, rhou, rhoE, m, S, L, v, gamma, CFL, t)
    N = length(rho);
    h = 1 / (N-1);
    y = 0:h:1;
    
    u = rhou./rho;
    E = rhoE./rho;
    e = E - u.^2/2;
    
    p = (gamma-1).*rho.*e;
    
    % boundary
    p(1) = p(2);
    rho(1) = rho(2);
    e(1) = e(2);
    u(1) = 0;
    
    p(N) = p(N-1);
    rho(N) = rho(N-1);
    e(N) = e(N-1);
    u(N) = v;
    
    E = e + u.^2/2;
    
    c_sound = max(sqrt(gamma*p./rho));
    vmax = max(c_sound,max(u));
    dt = CFL * h * L / vmax;
    
    % mass flux
    f1 = rho.*u;
    % momentum flux
    f2 = rho.*u.^2 + p;
    % energy flux
    f3 = rho.*u.*E + p.*u;
    
    % pseudo mass flux
    f11 = rho;
    % pseudo momentum flux
    f22 = rho.*u;
    % pseudo energy flux
    f33 = rho.*E;
    
    i = 2:N-1;
    im = 1:N-2;
    ip = 3:N;
    
    rho1 = rho;
    rhou1 = rhou;
    rhoE1 = rhoE;
    
    rho1(i) = rho(i) - dt / (L * h * 2) * (f1(ip)-f1(im)) - dt*v*y(i)/(L * h) .* (f11(i)-f11(im));
    rhou1(i) = rhou(i) - dt / (L * h * 2) * (f2(ip)-f2(im)) - dt*v*y(i)/(L * h) .* (f22(i)-f22(im));
    rhoE1(i) = rhoE(i) - dt / (L * h * 2) * (f3(ip)-f3(im)) - dt*v*y(i)/(L * h) .* (f33(i)-f33(im));
    
    rho1 = filter_neumann(rho1);
    rhou1 = filter_neumann(rhou1);
    rhoE1 = filter_neumann(rhoE1);
    
    L1 = L + v * dt;
    v1 = v + p(N) * S / m * dt;   
    
    t1 = t + dt;
    
    mass = L * h * sum(rho(i));
end