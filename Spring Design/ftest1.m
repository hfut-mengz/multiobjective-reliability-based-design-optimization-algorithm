function [f,g,rst] = ftest1(x,varargin)
% 10-bar 2D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = [0.51,0.6,0.51]';% lower-bound of design variables
    UB = [70.49,3,42.49]';% upper-bound of design variables
    
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x1 = round(x(1));
    x2 = x(2);
    d  = [0.009,0.0095,0.0104,0.0118,0.0128,0.0132,0.014,....
              0.015, 0.0162, 0.0173, 0.018, 0.020, 0.023, 0.025,...
              0.028, 0.032, 0.035, 0.041, 0.047, 0.054, 0.063,....
              0.072, 0.080, 0.092, 0.0105, 0.120, 0.135, 0.148,....
              0.162, 0.177, 0.192, 0.207, 0.225, 0.244, 0.263,....
              0.283, 0.307, 0.331, 0.362,0.394,0.4375,0.500];
   x3 = d(max(1,min(42,round(x(3))))); 
   %%%%%%%%%
   %% constants
        cf = (4.*x2./x3-1)./(4.*x2./x3-4)+0.615.*x3./x2;
        K  = (11.5.*10.^6.*x3.^4)./(8.*x1.*x2.^3);
        lf = 1000./K + 1.05.*(x1+2).*x3;
        sigp = 300./K;
   %% objective function
        f(1) = (pi.^2.*x2.*x3.^2.*(x1+2))./4;
        f(2) = (8000.*cf.*x2)./(pi.*x3.^3);
    %%%%%%constraint
    nr=3;
    ns=3;
    ng=8; %num of constraint
    betat=3;
    [fval,gnum,dgnum] = PMA(x,betat,nr,ns,ng);
    g(:,1)=fval;
        
    % Solving with Finit Element Method (FEM)
    rst=[];
    
end