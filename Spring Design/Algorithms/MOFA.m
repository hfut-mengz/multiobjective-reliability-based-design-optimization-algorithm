%% This demo shows how the multiobjective firefly algorithm (MOFA) works  %
%% The standard firefly algorithm (FA) was developed by Xin-She Yang and  %
%% its Matlab code is alrealy available online at Mathworks, etc.         %
%% This demo focuses on the MOFA with non-dominated sorting               %
%% ---------------------------------------------------------------------- % 

%% Programmed by Xin-She Yang in 2011 and 2012 @ Cambridge Univ and 
%% National Physics Laboratory, London. Updated and last modified 
%% by X S Yang in 2014, and then 2015. 
% References:
% (1) Xin-She Yang, Firefly algorithm, stochastic test functions and 
%     design optimizaiton, Int. Journal of Bio-Inspired Computation, 
%     vol. 2, no. 2, 78-84 (2010). 
% (2) Xin-She Yang, Multiobjective firely algorithm for continuous
%     optimizatoin, Engineering with Computers, vol. 29, no. 2, 
%     175--184 (2013).
% (3) Xin-She Yang, Nature-Inspired Optimization Algorithms, 
%     Elsevier Insight, (2014).  [Book]
% -----------------------------------------------------------------------

function rst=MOFA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
% if nargin<1,    
%   inp=[100 1000]; % Default parameters
% end
n=nsol;               % Population size (number of fireflies)
tMax=nloop;            % Maximum number of iterations
alpha=1.0;              % Randomness strength 0--1 (highly random)
beta0=1.0;              % Attractiveness constant
gamma=0.1;              % Absorption coefficient
theta=10^(-4/tMax);     % The parameter theta can be taken as 0.97 to 0.99    
                        % This is a randomness reduction factor for alpha
% G=zeros(nloop,nvar)';
% For the ZDT Function #3 with m=2 objectives
% m=3;             % Number of objectives
RnD=zeros(n,2);  % Initilize the rank and distance matrix
% Dimension of the search/independent variables
d=nvar;
Lb=a';   % Lower bounds/limits
Ub=b';   % Upper bounds/limits
%%number of objective
x0=Lb+(Ub-Lb).*rand(1,d);
[f0,G0]=feval(fun,x0');
f0=fpenal_1(f0,G0);
m=length(f0);
% Generating the initial locations of n fireflies
for i=1:n,
   Sol(i,:)=Lb+(Ub-Lb).*rand(1,d); 
   [f(i,1:m),G(i,:)] = feval(fun,Sol(i,:)');
   f(i,1:m)=fpenal_1(f(i,1:m),G(i,:));
end
% Store the fitness or objective values
f_new=f;
%% Sort the initialized population
x=[Sol f];  % combined into a single input
% Non-dominated sorting for the initila population
Sorted=solutions_sorting(x, m,d);
% Decompose into solutions, fitness, rank and distances
Sol=Sorted(:,1:d);  S_new=Sol;       % Record solutions 
f=Sorted(:,(d+1):(d+m));  f_new=f;   % Record objectives
RnD=Sorted(:,(d+m+1):end);           % Record ranks
   
for t=1:tMax,          %%%%% start the firely algorithm iterations %%%%%
   alpha=alpha*theta;  % Reduce alpha by a factor 0<theta<1
   scale=abs(Ub-Lb);   % Scale of the optimization problem
   Sol_old=Sol;        % Save the old population
   f_old=f;            % Save the old population objectives
% Two loops over all the n fireflies
for i=1:n,
    for j=i:n, 
      % Update moves and move to the brighter/more attractive
      % That is, all m objectives [i.e., f(,1:m)] should improve. 
      % For example, for m=2, this means that the logical 
      % condition (f(j,1)<=f(i,1) & f(j,2) <=f(i,2)) is true.     
      if (f(j,1:m)<=f(i,1:m)),     
         r=sqrt(sum((Sol(i,:)-Sol(j,:)).^2));
         beta=beta0*exp(-gamma*r.^2);     % Attractiveness
         steps=alpha.*(rand(1,d)-0.5).*scale;
      % The FA equation for updating position vectors 
      % That is, to move firefly i torwards firefly j
         Sol(i,:)=Sol(i,:)+beta*(Sol(j,:)-Sol(i,:))+steps;
         Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
      end
         [f(i,1:m),G(i,:)]=feval(fun,Sol(i,1:d)'); 
         f(i,1:m)=fpenal_1(f(i,1:m),G(i,:));
   end % end for j
end % end for i
   
   %% Evalute the fitness/function values of the new population
   for j=1:n,
        [f_new(j, 1:m),G(j, :)]=feval(fun,Sol(j,1:d)');
        f_new(j,1:m)=fpenal_1(f_new(j, 1:m),G(j,:));
        if (f_new(j,1:m) <= f(j,1:m)),   % if all improve
            f(j,1:m)=f_new(j,1:m);
        end
        % Update the current best (stored in the first row)
        if (f_new(j,1:m) <= f(1,1:m)), 
            Sol(1,1:d) = Sol(j,1:d); 
            f_new(1,:)=f_new(j,:);
        end
    end % end of for loop j
     
%% ! It's very important to combine both populations, otherwise,
%% the results may look odd and will be very inefficient. !     
%% The combined population consits of both the old and new solutions
%% So the total size of the combined population for sorting is 2*n
       X(1:n,:)=[Sol f_new];               % Combine new solutions
       X((n+1):(2*n),:)=[Sol_old f_old];   % Combine old solutions
       Sorted=solutions_sorting(X, m, d);
       %% Select n solutions among a combined population of 2*n solutions
       new_Sol=Select_pop(Sorted, m, d, n);
       % Decompose into solutions, fitness and ranking
       Sol=new_Sol(:,1:d);             % Sorted solutions
       f=new_Sol(:,(d+1):(d+m));       % Sorted objective values
       RnD=new_Sol(:,(d+m+1):end);     % Sorted ranks and distances
    
%   %% Running display at each 100 iterations
%    if ~mod(t,10), 
%      disp(strcat('Iterations t=',num2str(t))); 
%      plot(f(:, 1), f(:, 2),'ro','MarkerSize',3); 
%      axis([0 1 -0.8 1]);
%      xlabel('f_1'); ylabel('f_2');
%      drawnow;
%    end   
%     [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],Sol',f',G',narchive);
    rst.ppareto{t}=Sol';
    rst.fpareto{t}=f';
    rst.gpareto{t}=G';
    rst.timestamp=datetime('now');

end % End of t loop (up to tMax) and end of the main FA loop  

%% Make sure that new fireflies are within the bounds/limits
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective functions
function f = obj_funs(x, m)
% Zitzler-Deb-Thiele's funciton No 3 (ZDT function 3)
% m = # of objectives   % d = # of variables/dimensions
d=length(x);  % d=30 for ZDT 3
% First objective f1
f(1) = x(1);
g=1+9/29*sum(x(2:d));
h=1-sqrt(f(1)/g)-f(1)/g*sin(10*pi*f(1));
% Second objective f2
f(2) = g*h;

%%%%%%%%%%%%%%%%%% end of the definitions of obojectives %%%%%%%%%%%%%%%%%%

function new_Sol = Select_pop(firefly, m, ndim, npop)
% The input population to this part has twice (ntwice) of the needed 
% population size (npop). Thus, selection is done based on ranking and 
% crowding distances, calculated from the non-dominated sorting
ntwice= size(firefly,1);
% Ranking is stored in column Krank
Krank=m+ndim+1;
% Sort the population of size 2*npop according to their ranks
[~,Index] = sort(firefly(:,Krank)); 
sorted_firefly=firefly(Index,:);
% Get the maximum rank among the population
RankMax=max(firefly(:,Krank)); 

%% Main loop for selecting solutions based on ranks and crowding distances
K = 0;  % Initialization for the rank counter 
% Loop over all ranks in the population
for i =1:RankMax,  
    % Obtain the current rank i from sorted solutions
    RankSol = max(find(sorted_firefly(:, Krank) == i));
    % In the new solutions, there can be npop solutions to fill
    if RankSol<npop,
       new_Sol(K+1:RankSol,:)=sorted_firefly(K+1:RankSol,:);
    end 
    % If the population after addition is large than npop, re-arrangement
    % or selection is carried out
    if RankSol>=npop
        % Sort/Select the solutions with the current rank 
        candidate_firefly=sorted_firefly(K + 1 : RankSol, :);
        [~,tmp_Rank]=sort(candidate_firefly(:,Krank+1),'descend');
        % Fill the rest (npop-K) fireflies/solutions up to npop solutions 
        for j = 1:(npop-K), 
            new_Sol(K+j,:)=candidate_firefly(tmp_Rank(j),:);
        end
    end
    % Record and update the current rank after adding new solutions
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