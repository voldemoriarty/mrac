function [alpha, beta, x0] = linearize(y0)
    a = 1.65;
    b = 6.2;
    m = 0.12;
    g = 9.8;
    
    % compute x0 corresponding to the height given
    % use the nonlinear height to force relation
    x0 = a * ((y0 + b)^4) * m * g;
    
    % now compute alpha and beta, derivation in the 
    % notes
    alpha = 4*x0/(a*((y0+b)^5));
    beta  = 1/(a*((y0+b)^4)); 
end