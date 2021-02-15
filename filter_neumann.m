function [g] = filter_neumann(f)
    N = length(f);
    i = 2:N-1;
    im = 1:N-2;
    ip = 3:N;
    
    f1 = f;
    f1(1) = f1(2);
    f1(N) = f1(N-1);
    
    g = f1;
    g(i) = (f1(im) + 2 * f1(i) + f1(ip)) / 4;    
end

