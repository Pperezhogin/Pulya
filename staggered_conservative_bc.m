function [rho1,rhou1,rhoE1,L1,v1,t1,dt,u,p,e,c_sound,mass] = staggered_conservative_bc(rho, rhou, rhoE, m, S, L, v, gamma, CFL, t)
    % model variables are in the center of the cells
    N = length(rho);
    h = 1 / (N);
    y = h/2:h:1-h/2;
    
    u = rhou./rho;
    E = rhoE./rho;
    e = E - u.^2/2;
    
    p = (gamma-1).*rho.*e;
   
    c_sound = sqrt(gamma*p./rho);
    vmax = max(max(c_sound),max(u));
    dt = CFL * h * L / vmax;
    
    % mass flux
    f1 = rhou + v*rho.*y;
    % momentum flux
    f2 = rhou.*u + p + v*rhou.*y;
    % energy flux
    f3 = rhou.*E + p.*u + v*rhoE.*y;
    
    % in y in range 0:h:1
    f11 = zeros(1,N+1);
    f22 = zeros(1,N+1);
    f33 = zeros(1,N+1);
    im = 1:N-1;
    ip = 2:N;
    i = 2:N;
    f11(i) = (f1(im) + f1(ip)) / 2;
    f22(i) = (f2(im) + f2(ip)) / 2;
    f33(i) = (f3(im) + f3(ip)) / 2;
    
    f11(1) = 0;
    f22(1) = p(1)*0;
    f33(1) = 0;
    f11(N+1) =  2 * v * h * sum(rho);
    f22(N+1) =  2 * v * h * sum(rhou) + p(N)*0;
    f33(N+1) =  2 * v * h * sum(rhoE) + p(N) * v*0;
    
    ip = 2:N+1;
    im = 1:N;
    rho1 = rho - dt / (L * h) * (f11(ip)-f11(im)) + dt*v/L*rho;
    rhou1 = rhou - dt / (L * h) * (f22(ip)-f22(im)) + dt*v/L*rhou;
    rhoE1 = rhoE - dt / (L * h) * (f33(ip)-f33(im)) + dt*v/L*rhoE;
    
    rho1 = filter_neumann(rho1);
    rhou1 = filter_neumann(rhou1);
    rhoE1 = filter_neumann(rhoE1);
    
    L1 = L + v * dt;
    v1 = v + p(N) * S / m * dt;   
    
    t1 = t + dt;
    
    mass = L * h * sum(rho);
end