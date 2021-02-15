function [rho1,rhou1,rhoE1,L1,v1,t1,dt,u,p,e,c_sound,mass] = collocated_conservative_bc(rho, rhou, rhoE, m, S, L, v, gamma, CFL, t)
    N = length(rho);
    h = 1 / (N-1);
    y = 0:h:1;
    i = 2:N-1;
    im = 1:N-2;
    ip = 3:N;
    
    u = rhou./rho;
    E = rhoE./rho;
    e = E - u.^2/2;
    
    p = (gamma-1).*rho.*e;
    
    % boundary
    p(1) = p(2);    
    p(N) = p(N-1);
    u(1) = 0;
    u(N) = v;
    
    c_sound = sqrt(gamma*p./rho);
    vmax = max(max(c_sound(i)),max(u(i)));
    dt = CFL * h * L / vmax;
    
    % mass flux
    f1 = rhou;
    % momentum flux
    f2 = rhou.*u + p;
    % energy flux
    f3 = rhou.*E + p.*u;
    
    % B.C.
    f1(1) = 0;
    f1(N) = v * h * sum(rho(i));
    f2(1) = p(1);
    f2(N) = v * h * sum(rhou(i)) + p(N);
    f3(1) = 0;
    f3(N) = v * h * sum(rhoE(i)) + p(N) * v;
    
    
    % pseudo mass flux
    f11 = v*rho.*y;
    % pseudo momentum flux
    f22 = v*rhou.*y;
    % pseudo energy flux
    f33 = v*rhoE.*y;
    % B.C.
    f11(1) = 0;
    f22(1) = 0;
    f33(1) = 0;
    f11(N) = v * h * sum(rho(i));
    f22(N) = v * h * sum(rhou(i));
    f33(N) = v * h * sum(rhoE(i));    
    
    rho1 = rho;
    rhou1 = rhou;
    rhoE1 = rhoE;
    
    rho1(i) = rho(i) - dt / (L * h * 2) * (f1(ip)-f1(im)) - dt/(L * h * 2) * (f11(ip)-f11(im)) + dt*v/L*rho(i);
    rhou1(i) = rhou(i) - dt / (L * h * 2) * (f2(ip)-f2(im)) - dt/(L * h * 2) * (f22(ip)-f22(im)) + dt*v/L*rhou(i);
    rhoE1(i) = rhoE(i) - dt / (L * h * 2) * (f3(ip)-f3(im)) - dt/(L * h * 2) * (f33(ip)-f33(im)) + dt*v/L*rhoE(i);
    
    rho1 = filter_neumann(rho1);
    rhou1 = filter_neumann(rhou1);
    rhoE1 = filter_neumann(rhoE1);
    
    L1 = L + v * dt;
    v1 = v + p(N) * S / m * dt;   
    
    t1 = t + dt;
    
    mass = L * h * sum(rho(i));
end