clear all;close all;warning off;clc;
addpath([pwd '\Algorithms']);
%%%%%%%%%%%%%%
%   INPUT    %
%%%%%%%%%%%%%%

nloop=100;      % Number of generations
nsol=100;       % Population size


narchive=nsol;  % Pareto archive size
Nrun=20;        % Number of runs per algorithm per problem
%%%%%%%%%%%%%%
% algo={'MOALO'; 'MODA';'MOGOA';'MOGWO'; 'MOMVO';
%         'MOWCA'; 'MSSA';'SHAMODE';'SHAMODE_WO';
%         'NSGA_II';'RPBILDE';'DEMO';'MOEA_D';'UPSEMOA'};
algo={'NSGA_II';'DEMO';'MOEA_D';'MOPSO';'MOFPA';'MOBA_NS';...
    'MOGWO'; 'MOMVO';'MOWCA';'SHAMODE';'SHAMODE_WO';...
    'MSSA';'RPBILDE';'UPSEMOA';'NSMFO'};                                            
    
fobj={  'ftest1'};    
nvar=[4];                % no. of design variables of the test problems

for k=1:numel(fobj) %%%% benchmark
    nbit=1000;
    funj=char(fobj(k,:));
    nvari=nvar(k);
    % Lower and upper bound of all design variables are set as 0 and 1
    % respectively. Design variables will be modified to real values inside
    % objective function files (fxxbarSI.m).
    a=zeros(nvari,1);
    b=ones(nvari,1);
       
    for i=1:numel(algo)  %%%% algorithm
        filename=['rst_' num2str(k,'%03.f') '_' num2str(i,'%03.f') '.mat'];
        for j=1:Nrun     %%%%   optimization run
            RR=[k i j]
            rst(j)=feval(char(algo(i,:)),funj,filename,nloop,nsol,nvari,nbit,narchive,a,b);
        end
        save(filename,'rst','-v7.3');
    end
end   