%% Multiobjective Bat Algorithm (MOBA) with non-dominated sorting (NS).   %
%% MOBA was developed by Xin-She Yang in 2011.                            %
%% The original verison of MOBA without non-dominated sorting was         %
%% implemented by Xin-She Yang in Nov 2010, then finalized in Jan 2011.   %
%% This version of MOBA has been implemented as a combination of MOBA     %
%% with non-dominated sorting. The implementation was based on the code   %
%% by M. Jamil in July 2015 and earlier MOBA codes by X S Yang in 2011.   %
%% Updated and last modified by X. S. Yang in Sept 2015                   %                

%% Some Relevant References: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Xin-She Yang, Bat algorithm for multi-objective optimisaiton,       % 
% Int. Journal of Bio-Inspired Computation, vol. 3, no.4, 267-274 (2011). %
% (2) M. Jamil, H.J. Zepernick, X.S. Yang, Synthesizing cross-ambiguity   %
% functions using the improved bat algorithm, Recent Advanceds in Swarm   % 
% Intelligence and Evolutionary Computation (Ed. X.S. Yang), Springer,    %
% pp.179-202 (2015). [book chapter.]                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rst = MOBA_NS(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin<1,      % Check inputs (otherwise, default values)
% inp=[100 1000];   % [Pop size, #Iteration] e.g., inp=[100 1000]
% end
% G=zeros(1,nsol);
npop   = nsol;  %% Population size
Gen = nloop;     %% Number of Iterations
m=3;              %% Number of objectives
ndim=nvar;          %% ndim=dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter settings  (initial values)                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AL      = 0.9; 
r       = 0.9; 
alpha   = 0.9; 
Gamma   = 0.9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Initialize the increment components for bats/locations
del_bats = zeros(npop,ndim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_bats = zeros(npop,ndim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lb=a';     % Lower bounds
Ub=b';      % Upper bounds
minf=0; maxf=1;       % Frequnecy range

x0=Lb+(Ub-Lb).*rand(1,ndim);
[f0,G0]=feval(fun,x0');

m=length(f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the population
for i=1:npop
   x(i,:)=Lb+(Ub-Lb).*rand(1,ndim); 
   [f(i,1:m),G(:,i)] = feval(fun,x(i,:)');
   f(i,1:m)=fpenal_1(f(i,1:m),G(:,i));
   
end
%% Combining x and f into a large matrix for easy storage
%% x is an npop x ndim matrix, while f is an (npop x m) matrix.
%% So the combined matrix has a size of npop x (ndim+m),
%% with the last m columns being the m objective values
bats=[x f]; 

%% Sort the initialized population
% Sort the population using nondomination sorting. This returns two columns
% for each individual solution, and these 2 columns correspond to the rank 
% and the crowding distance (based on their positions in the front). 
bats = solutions_sorting(bats, m, ndim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : Gen, 
    %%%%%  Update loudness and pulse emission rate %%%%%%%%%%%%%%%%%%%%%%%%
    AL = (AL ^ i)*alpha;
    r = r*(1 - Gamma ^ i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:npop, 
        %% Varying frequency
        freq(j) = minf + (maxf - minf) * rand;
        %% Updating velocity
        %% Notes: 
        %% The best bat based on crowding dist is stored in Bats(1,1:ndim)
        best=bats(1,1:ndim);
        del_bats(j,:)=del_bats(j,:)+(bats(j,1:ndim)-best).*freq(j);
        %% Update locations of bats (see Paper by X.S. Yang, 2010)
        new_bats(j,1:ndim) = bats(j,1:ndim) + del_bats(j,1:ndim);
        if (rand > r)
            %% Randomization and generation of new bats
            epsilon=randn(1,ndim); 
            dS = epsilon.*(bats(j,1:ndim) - bats(1,(1:ndim)));
            new_bats(j,1:ndim) = bats(j,1:ndim) + dS;
        end
        %% Check if new solutions are within limits
        new_bats(j,1:ndim) = findlimits(new_bats(j,1:ndim),Lb, Ub);
        %% Evalute the fitness/function values of the new population
        [new_bats(j, ndim+1:m+ndim),G(:,j)] = feval(fun,new_bats(j,1:ndim)');%obj_funs(new_bats(j,1:ndim),m);
        new_bats(j, ndim+1:m+ndim)=fpenal_1(new_bats(j, ndim+1:m+ndim),G(:,j));
        if ((new_bats(j,ndim+1:m+ndim)<=bats(j,ndim+1:m+ndim))|(rand<AL))
            bats(j,ndim:ndim+1)  = new_bats(j,ndim:ndim+1);
            bats(j,ndim+1:m+ndim)= new_bats(j,ndim+1:m+ndim);
        end
        % Update the current best bat (stored in the first row)
        if new_bats(j,ndim+1:ndim+m) <= bats(1,(ndim+1:ndim+m)) 
            bats(1,1:(ndim+m)) = new_bats(j,1:(ndim+m)); 
        end
    end
    
%% The combined population consits of both old bats and new bats
%% So the total population for sorting is 2*npop 
%% ! Very important to combine old and new bats !
   Sort_bats(1:npop,:) = bats;
   Sort_bats((npop + 1):(2*npop), 1:m+ndim) = new_bats;
%% Non-dominated sorting process (a separate function/subroutine)
   Sorted_bats = solutions_sorting(Sort_bats, m, ndim);
 
%% Select npop solutions among a combined population of 2*npop solutions  
    bats = cleanup_batspop(Sorted_bats, m, ndim, npop);  

%% Running display at each 100 iterations
%    if ~mod(i,10), 
%       disp(strcat('Iterations t=',num2str(i))); 
%       plot(bats(:, ndim+1), bats(:, ndim+2),'ro','MarkerSize',5); 
%       xlabel('f_1'); ylabel('f_2');
%       drawnow;
%    end
   
%    [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],bats(:,1:ndim)',bats(:,ndim+1:ndim+m)',G,narchive);
    rst.ppareto{i}=bats(:,1:ndim)';
    rst.fpareto{i}=bats(:,ndim+1:ndim+m)';
    rst.gpareto{i}=G;
    rst.timestamp=datetime('now');
end

%% Displaying the final Pareto Front as a graph
% figure(2)    
% plot(bats(:,ndim+1), bats(:,ndim+2),'rd');
% title(strcat(num2str(npop),' Points on the Pareto Front'));
% xlabel('Objective (f_1)'); ylabel('Objective (f_2)');


    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure all the bats are within the limits
function [ns]=findlimits(ns,Lb,Ub)
  % Apply the lower bound
  ns_tmp=ns;
  I=ns_tmp < Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  ns=ns_tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective functions
% function f = obj_funs(x, m)
% % Zitzler-Deb-Thiele's funciton No 3 (ZDT function 3)
% % m    = number of objectives
% % ndim = number of variables/dimensions
% ndim=length(x);  % ndim=30 for ZDT 3
% % First objective f1
% f(1) = x(1);
% g=1+9/29*sum(x(2:ndim));
% h=1-sqrt(f(1)/g)-f(1)/g*sin(10*pi*f(1));
% % Second objective f2
% f(2) = g*h;
%%%%%%%%%%%%%%%%%% end of the definitions of obojectives %%%%%%%%%%%%%%%%%%

%% Clean up the populations (both old and new) to give a new population
% This cleanup here is similar to the Non-dominated Sorting Genetic
% Algorithm (NSGA-II) by K. Deb et al. (2002), which can be applied to 
% any cleanup of 2*npop solutions to form a set of npop solutions.
function new_bats = cleanup_batspop(bats, m, ndim, npop)
% The input population to this part has twice (ntwice) of the needed 
% population size (npop). Thus, selection is done based on ranking and 
% crowding distances, calculated from the non-dominated sorting
ntwice= size(bats,1);
% Ranking is stored in column Krank
Krank=m+ndim+1;
% Sort the population of size 2*npop according to their ranks
[~,Index] = sort(bats(:,Krank)); sorted_bats=bats(Index,:);
% Get the maximum rank among the population
RankMax=max(bats(:,Krank)); 

%% Main loop for selecting solutions based on ranks and crowding distances
K = 0;  % Initialization for the rank counter 
% Loop over all ranks in the population
for i =1:RankMax,  
    % Obtain the current rank i from sorted solutions
    RankSol = max(find(sorted_bats(:, Krank) == i));
    % In the new bats/solutions, there can be npop solutions to fill
    if RankSol<npop,
       new_bats(K+1:RankSol,:)=sorted_bats(K+1:RankSol,:);
    end 
    % If the population after addition is large than npop, re-arrangement
    % or selection is carried out
    if RankSol>=npop
        % Sort/Select the solutions with the current rank 
        candidate_bats = sorted_bats(K + 1 : RankSol, :);
        [~,tmp_Rank]=sort(candidate_bats(:,Krank+1),'descend');
        % Fill the rest (npop-K) bats/solutions up to npop solutions 
        for j = 1:(npop-K), 
            new_bats(K+j,:)=candidate_bats(tmp_Rank(j),:);
        end
    end
    % Record and update the current rank after adding new bats 
    K = RankSol;
end

function sorted_x = solutions_sorting(x, m, ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs and outputs are the extended solutions x with a dimension of    %
%% npop by (ndim+m+2). The objective values are already included in x.    %
% More specifically, the first ndim columns are the actual solutions or 
% variable values (1:ndim), followed by the m columns of objective values. 
% Then, the next column (i.e.,ndim+m+1) corresponds to the ranks, whereas 
% the final column (i.e., ndim+m+2) records the crowd distances. 
% ----------------------------------------------------------------------- %

%% Get the parameters from the inputs such as the size of input population
npop=size(x,1);    % Population size
frontRank=1;       % Pareto frontRank (counter) initialization
Rcol=ndim+m+1;     % Store the ranks in the column Rcol=ndim+m+1
% Define the Parato Front as a class (PF) and initilization of xSol
PF(frontRank).R=[];   xSol=[];
%% The main non-dominated sorting starts here             %%%%%%%%%%%%%%%%% 
for i = 1:npop, 
    % Set the number (initially, 0) of solutions dominating this solution
    xSol(i).n=0;
    % Find all the solutions (that dominated by this solution)
    xSol(i).q=[];
    % Sorting into 3 categories: better (minimization), equal & otherwise
    for j=1:npop,
        % Definte 3 counters for 3 categories
        ns_categ_1=0; ns_categ_2=0; ns_categ_3=0;
        for k=1:m,  % for all m objectives
            % Update the counters for 3 different categories
            if (x(i,ndim+k) < x(j,ndim+k)),      % better/non-dominated
                ns_categ_1=ns_categ_1+1;
            elseif (x(i,ndim+k)==x(j,ndim+k)),   % equal
                ns_categ_2=ns_categ_2+1;
            else                                 % dominated
                ns_categ_3=ns_categ_3+1;
            end
        end % end of k
        % Update the solutions in their class
        if ns_categ_1==0 && ns_categ_2 ~= m
            xSol(i).n=xSol(i).n+1;
        elseif ns_categ_3 == 0 && ns_categ_2 ~= m
            xSol(i).q=[xSol(i).q j];
        end
    end % end of j   
    %% Record/Udpate the Pareto Front
    if xSol(i).n==0, 
        x(i,Rcol)=1;   % Update the Rank #1 (i.e., the Pareto Front)
        PF(frontRank).R = [PF(frontRank).R i];
    end
end % end of i=1:npop (The first round full rank-sorting process)

% Update the rest frontRanks (close, but not on the Pareto Front)
while ~isempty(PF(frontRank).R),
    nonPF=[];    % Intialization the set
    N=length(PF(frontRank).R);
for i=1 :N, 
   % Get the solution/list 
   Sol_tmp_q=xSol(PF(frontRank).R(i)).q; 
   % If not empty, update 
   if ~isempty(xSol(Sol_tmp_q))
       for j = 1:length(Sol_tmp_q),
         % Get the solutions dominated by the current solution    
          Sol_tmp_qj=xSol(PF(frontRank).R(i)).q(j);   
          xSol(Sol_tmp_qj).n=xSol(Sol_tmp_qj).n-1;
          if xSol(Sol_tmp_qj).n==0
             x(Sol_tmp_qj, Rcol)=frontRank + 1;
             nonPF = [nonPF Sol_tmp_qj];
          end
       end % end of j
   end
end  % end of i
   frontRank=frontRank+1;
   PF(frontRank).R=nonPF;
end % end of PF(frontRank)

% Now carry out the sorting of ranks and then update 
[~,frontRanks_Index]=sort(x(:, Rcol));
Sorted_frontRank=x(frontRanks_Index,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate the crowding distances for each solution for each frontRank   %
% That is, all the non-domonated solutions on the Pareto Front.  %%%%%%%%%%
Qi=0;      % Initialize a counter
for frontRank=1:(length(PF)-1), 
    % Define/initialize a generalized distance matrix 
    dc = [];    past_Q=Qi+1;
    for i=1:length(PF(frontRank).R),
        dc(i,:)=Sorted_frontRank(Qi+i,:);
    end
    Qi=Qi+i;
    % Solutions are sorted according to their fitness/objective values
    fobj_sorted=[];
    for i=1:m, 
        [~, f_Rank]=sort(dc(:,ndim+i));
        fobj_sorted=dc(f_Rank,:);
        % Find the max and min of the fobj values   
        fobj_max=fobj_sorted(length(f_Rank), ndim+i);
        fobj_min=fobj_sorted(1, ndim+i);
        % Calculate the range of the fobj
        f_range=fobj_max-fobj_min;
        % If the solution is at the end/edge, set its distance as infinity
        dc(f_Rank(length(f_Rank)), Rcol+i)=Inf;
        dc(f_Rank(1), Rcol+i) = Inf;
        for j=2:length(f_Rank)-1, 
            fobj2=fobj_sorted(j+1,ndim + i);
            fobj1=fobj_sorted(j-1,ndim + i);  
            % Check the range or special cases
            if (f_range==0),
                dc(f_Rank(j), Rcol+i)=Inf;
            else
            % Scale the range for distance normalization     
            dc(f_Rank(j),Rcol+i)=(fobj2-fobj1)/f_range;
            end
        end % end of j
    end % end of i
    
    % Calculate and update the crowding distances on the Pareto Front
    dist = []; dist(:,1)=zeros(length(PF(frontRank).R),1);
    for i=1:m, 
        dist(:,1)=dist(:,1)+dc(:, Rcol+i);
    end
    % Store the crowding distrance (dc) in the column of Rcol+1=ndim+m+2
    dc(:, Rcol+1)=dist;  dc=dc(:,1:Rcol+1);
    % Update for the output
    xy(past_Q:Qi,:)=dc;  
end  % end of all ranks search/update
sorted_x=xy();    % Output the sorted solutions

function [pareto1,fpareto1,gpareto1,A]=pbil_selection0(x1,f1,g1,pareto,fpareto,gpareto,narchive)

    %
    % GA Elite strategy
    % keep 1 elite from the old generation and another from the
    % new generation
    x=[pareto x1];
    f=[fpareto f1];
    g=[gpareto g1];

    [m0,n0]=size(fpareto);
    [m1,n1]=size(x);
    
    for i=1:n1
        xi=x(:,i);
        fi=f(:,i);
        gi=g(:,i);
        A(i,i)=0;
        for j=(i+1):n1
            xj=x(:,j);
            fj=f(:,j);
            gj=g(:,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
            A(i,j)=p_count1;
            A(j,i)=p_count2;
            %%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B=sum(A,1);
    Indm=[];
    for i=1:n1
        if B(i)==0
            Indm=[Indm i];
        end
    end
    nndm=length(Indm);

    pareto1=x(:,Indm);
    fpareto1=f(:,Indm);
    gpareto1=g(:,Indm);
    A=A(Indm,Indm);


function [p1,p2]=fdominated(f1,g1,f2,g2)
    n=length(f1);
    mg1=max(g1);
    mg2=max(g2);

    icount11=0;
    icount12=0;
    icount21=0;
    icount22=0;

    if mg1<=0 && mg2<=0
        for i=1:n
            if f1(i) <= f2(i)
                icount11=icount11+1;
            end
            if f1(i) < f2(i)
                icount12=icount12+1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if f2(i) <= f1(i)
                icount21=icount21+1;
            end
            if f2(i) < f1(i)
                icount22=icount22+1;
            end
        end
        if icount11 == n && icount12 > 0
            p1=1;
        else
            p1=0;
        end
        if icount21 == n && icount22 > 0
            p2=1;
        else
            p2=0;
        end
    elseif mg1 <=0 && mg2 > 0
        p1=1;p2=0;
    elseif mg2 <=0 && mg1 > 0
        p1=0;p2=1;
    else
        if mg1 <= mg2
            p1=1;p2=0;
        else
            p1=0;p2=1;
        end
    end

function [fp] = fpenal_1(f,g)
    % penalty function
    if max(g)>0
        fp=f+100*max(g);
    else
        fp=f;
    end

