function [f,g,rst] = ftest1(x,varargin)
% 10-bar 2D truss optimization problem
% x = design variables (input - 0<=x<=1)
    
    % Decoding design variables
    LB = [10 10 0.9 0.9]';% lower-bound of design variables
    UB = [80 50 5 5]';% upper-bound of design variables
    x = LB+(UB-LB).*x; % truncate x to lower and upper-bound
    x1 = x(1);
	x2 = x(2);
	x3 = x(3);
	x4 = x(4);
   %%%%%%%%%
     P = 600;
     L = 200;
     E = 2e4;
	%% objectives
	f(1) = 2 .* x2 .* x4 + x3 .* (x1- 2.*x4);
	f(2) = P.*L.^3./(48.*E.*(x3 .*((x1 - 2.*x4).^3)+2.*x2.*x4.*(4.*x4.*x4+3.*x1.*(x1-2.*x4)))./12);
    %%%%%%constraint
    nr=4;
    ns=4;
    ng=1; %num of constraint
    betat=3;
    [fval,gnum,dgnum] = PMA(x,betat,nr,ns,ng);
    g(:,1)=fval;
        
    % Solving with Finit Element Method (FEM)
    rst=[];
    
end