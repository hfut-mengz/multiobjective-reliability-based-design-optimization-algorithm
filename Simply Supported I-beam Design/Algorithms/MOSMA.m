% M. Premkumar, P. Jangir, R. Sowmya, H. H. Alhelou, A. A. Heidari and H. Chen, 
% "MOSMA: Multi-objective Slime Mould Algorithm Based on Elitist Non-dominated Sorting," 
% in IEEE Access, doi: 10.1109/ACCESS.2020.3047936.
function rst = MOSMA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)   
dim=nvar;
N=nsol;
Max_iter=nloop;
X = zeros(nsol,nvar);
Sol = zeros(nsol,nvar);
weight = ones(nsol,nvar);%fitness weight of each slime mold
lb=a';
ub=b';
M=numel(feval(fun,lb)');
G=zeros(1,nsol);
%% Initialize the population
for i=1:nsol
   x(i,:)=lb+(ub-lb).*rand(1,dim); 
   f(i,1:M) = feval(fun,x(i,:));
end
new_Sol=[x f]; 
new_Sol = solutions_sorting(new_Sol, M, dim);
for i = 1 : nloop 
[SmellOrder,SmellIndex] = sort(Sol);  
worstFitness = SmellOrder(N);
bestFitness = SmellOrder(1);
S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
        for k=1:nsol
            if k<=(nsol/2)  
                weight(SmellIndex(k),:) = 1+rand()*log10((bestFitness-SmellOrder(k))/(S)+1);
            else
                weight(SmellIndex(k),:) = 1-rand()*log10((bestFitness-SmellOrder(k))/(S)+1);
            end
        end     
a = atanh(-(i/nloop)+1);   
b = 1-i/Max_iter;
    for j=1:nsol 
        best=(new_Sol(j,1:dim) - new_Sol(1,(1:dim)));
        if rand<0.03    
            X(j,:) = (ub-lb).*rand+lb;
        else
            p =tanh(abs(f(j)-best));  
            vb = unifrnd(-a,a,1,dim); 
            vc = unifrnd(-b,b,1,dim);        
                r = rand();
                A = randi([1,N]);  
                B = randi([1,N]);
                if r<p    
                    X(j,:) = best+ vb.*(weight(j,:).*X(A,:)-X(B,:));
                else
                    X(j,:) = best+ vc.*(weight(j,:).*X(A,:)-X(B,:));
                end  
        end
        Sol(j,1:dim) = X(j,1:dim);       
        Flag4ub=Sol(j,1:dim)>ub;
        Flag4lb=Sol(j,1:dim)<lb;
        Sol(j,1:dim)=(Sol(j,1:dim).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        %% Evalute the fitness/function values of the new population
        Sol(j, dim+1:M+dim) = feval(fun,Sol(j,1:dim));
        if Sol(j,dim+1:dim+M) <= new_Sol(1,(dim+1:dim+M)) 
           new_Sol(1,1:(dim+M)) = Sol(j,1:(dim+M));  
        end
    end    
%% ! Very important to combine old and new bats !
   Sort_bats(1:N,:) = new_Sol;
   Sort_bats((N + 1):(2*N), 1:M+dim) = Sol;
%% Non-dominated sorting process (a separate function/subroutine)
   Sorted_bats = solutions_sorting(Sort_bats, M, dim); 
%% Select npop solutions among a combined population of 2*npop solutions  
    new_Sol = cleanup_batspop(Sorted_bats, M, dim, N);  
    if rem(i, 10) == 0
    fprintf('Generation: %d\n', i);        
    end
    
    
    [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],new_Sol(:,1:dim)',new_Sol(:,dim+1:dim+M)',G,narchive);
    rst.ppareto{i}=ppareto;
    rst.fpareto{i}=fpareto;
    rst.gpareto{i}=gpareto;
    rst.timestamp=datetime('now');
end
f=new_Sol;
Obtained_Pareto= f(:,dim+1:dim+M); % extract data to plot
Obtained_Pareto=sortrows(Obtained_Pareto,3);
True_Pareto=load('ZDT3.txt');
M_IGD=IGD(Obtained_Pareto,True_Pareto);
M_GD=GD(Obtained_Pareto,True_Pareto);
M_HV=HV(Obtained_Pareto,True_Pareto);
M_Spacing=Spacing(Obtained_Pareto,True_Pareto);
M_Spread=Spread(Obtained_Pareto,True_Pareto);
M_DeltaP=DeltaP(Obtained_Pareto,True_Pareto);
display(['The IGD Metric obtained by MOSMA is     : ', num2str(M_IGD)]);
display(['The GD Metric obtained by MOSMA is      : ', num2str(M_GD)]);
display(['The HV Metric obtained by MOSMA is      : ', num2str(M_HV)]);
display(['The Spacing Metric obtained by MOSMA is : ', num2str(M_Spacing)]);
display(['The Spread Metric obtained by MOSMA is  : ', num2str(M_Spread)]);
display(['The DeltaP Metric obtained by MOSMA is  : ', num2str(M_DeltaP)]);
plotctrol(f,True_Pareto);

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

function Score = DeltaP(PopObj,PF)
% <metric> <min>
% Averaged Hausdorff distance

%------------------------------- Reference --------------------------------
% O. Schutze, X. Esquivel, A. Lara, and C. A. Coello Coello, Using the
% averaged Hausdorff distance as a performance measure in evolutionary
% multiobjective optimization, IEEE Transactions on Evolutionary
% Computation, 2012, 16(4): 504-522.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group.

    IGDp  = mean(min(pdist2(PF,PopObj),[],2));
    GDp   = mean(min(pdist2(PopObj,PF),[],2));
    Score = max(IGDp,GDp);

function Score = GD(PopObj,PF)
% <metric> <min>
% Generational distance

%------------------------------- Reference --------------------------------
% D. A. Van Veldhuizen, Multiobjective evolutionary algorithms:
% Classifications, analyses, and new innovations, Ph.D. thesis, Department
% of Electrical and Computer Engineering, Graduate School of Engineering,
% Air Force Institute of Technology, Wright Patterson Air Force Base, 1999.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. 
    Distance = min(pdist2(PopObj,PF),[],2);
    Score    = norm(Distance) / length(Distance);
    
    
 function [Score,PopObj] = HV(PopObj,PF)
% <metric> <max>
% Hypervolume

%------------------------------- Reference --------------------------------
% E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
% comparative case study and the strength Pareto approach, IEEE
% Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. 

    % Normalize the population according to the reference point set
    [N,M]  = size(PopObj);
    fmin   = min(min(PopObj,[],1),zeros(1,M));
    fmax   = max(PF,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat((fmax-fmin)*1.1,N,1);
    PopObj(any(PopObj>1,2),:) = [];
    % The reference point is set to (1,1,...)
    RefPoint = ones(1,M);
    if isempty(PopObj)
        Score = 0;
    elseif M < 4
        % Calculate the exact HV value
        pl = sortrows(PopObj);
        S  = {1,pl};
        for k = 1 : M-1
            S_ = {};
            for i = 1 : size(S,1)
                Stemp = Slice(cell2mat(S(i,2)),k,RefPoint);
                for j = 1 : size(Stemp,1)
                    temp(1) = {cell2mat(Stemp(j,1))*cell2mat(S(i,1))};
                    temp(2) = Stemp(j,2);
                    S_      = Add(temp,S_);
                end
            end
            S = S_;
        end
        Score = 0;
        for i = 1 : size(S,1)
            p     = Head(cell2mat(S(i,2)));
            Score = Score + cell2mat(S(i,1))*abs(p(M)-RefPoint(M));
        end
    else
        % Estimate the HV value by Monte Carlo estimation
        SampleNum = 1000000;
        MaxValue  = RefPoint;
        MinValue  = min(PopObj,[],1);
        Samples   = unifrnd(repmat(MinValue,SampleNum,1),repmat(MaxValue,SampleNum,1));
        if gpuDeviceCount > 0
            % GPU acceleration
            Samples = gpuArray(single(Samples));
            PopObj  = gpuArray(single(PopObj));
        end
        for i = 1 : size(PopObj,1)
            drawnow();
            domi = true(size(Samples,1),1);
            m    = 1;
            while m <= M && any(domi)
                domi = domi & PopObj(i,m) <= Samples(:,m);
                m    = m + 1;
            end
            Samples(domi,:) = [];
        end
        Score = prod(MaxValue-MinValue)*(1-size(Samples,1)/SampleNum);
    end


function S = Slice(pl,k,RefPoint)
    p  = Head(pl);
    pl = Tail(pl);
    ql = [];
    S  = {};
    while ~isempty(pl)
        ql  = Insert(p,k+1,ql);
        p_  = Head(pl);
        cell_(1,1) = {abs(p(k)-p_(k))};
        cell_(1,2) = {ql};
        S   = Add(cell_,S);
        p   = p_;
        pl  = Tail(pl);
    end
    ql = Insert(p,k+1,ql);
    cell_(1,1) = {abs(p(k)-RefPoint(k))};
    cell_(1,2) = {ql};
    S  = Add(cell_,S);


function ql = Insert(p,k,pl)
    flag1 = 0;
    flag2 = 0;
    ql    = [];
    hp    = Head(pl);
    while ~isempty(pl) && hp(k) < p(k)
        ql = [ql;hp];
        pl = Tail(pl);
        hp = Head(pl);
    end
    ql = [ql;p];
    m  = length(p);
    while ~isempty(pl)
        q = Head(pl);
        for i = k : m
            if p(i) < q(i)
                flag1 = 1;
            else
                if p(i) > q(i)
                    flag2 = 1;
                end
            end
        end
        if ~(flag1 == 1 && flag2 == 0)
            ql = [ql;Head(pl)];
        end
        pl = Tail(pl);
    end  


function p = Head(pl)
    if isempty(pl)
        p = [];
    else
        p = pl(1,:);
    end


function ql = Tail(pl)
    if size(pl,1) < 2
        ql = [];
    else
        ql = pl(2:end,:);
    end


function S_ = Add(cell_,S)
    n = size(S,1);
    m = 0;
    for k = 1 : n
        if isequal(cell_(1,2),S(k,2))
            S(k,1) = {cell2mat(S(k,1))+cell2mat(cell_(1,1))};
            m = 1;
            break;
        end
    end
    if m == 0
        S(n+1,:) = cell_(1,:);
    end
    S_ = S;     

    
    function Score = IGD(PopObj,PF)
% <metric> <min>
% Inverted generational distance

%------------------------------- Reference --------------------------------
% C. A. Coello Coello and N. C. Cortes, Solving multiobjective optimization
% problems using an artificial immune system, Genetic Programming and
% Evolvable Machines, 2005, 6(2): 163-190.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group.

    Distance = min(pdist2(PF,PopObj),[],2);
    Score    = mean(Distance);

    
    
function Score = Spacing(PopObj,PF)
% <metric> <min>
% Spacing
%------------------------------- Reference --------------------------------
% J. R. Schott, Fault tolerant design using single and multicriteria
% genetic algorithm optimization, Master's thesis, Department of
% Aeronautics and Astronautics, Massachusetts Institute of Technology,
% 1995.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. 
    Distance = pdist2(PopObj,PopObj,'cityblock');
    Distance(logical(eye(size(Distance,1)))) = inf;
    Score    = std(min(Distance,[],2));
    
    
function Score = Spread(PopObj,PF)
% <metric> <min>
% Spread

%------------------------------- Reference --------------------------------
% Y. Wang, L. Wu, and X. Yuan, Multi-objective self-adaptive differential
% evolution with elitist archive and crowding entropy-based diversity
% measure, Soft Computing, 2010, 14(3): 193-209.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group.

    Dis1  = pdist2(PopObj,PopObj);
    Dis1(logical(eye(size(Dis1,1)))) = inf;
    [~,E] = max(PF,[],1);PF(E,:)
    Dis2  = pdist2(PF(E,:),PopObj);
    d1    = sum(min(Dis2,[],2));
    d2    = mean(min(Dis1,[],2));
    Score = (d1+sum(abs(min(Dis1,[],2)-d2))) / (d1+(size(PopObj,1)-size(PopObj,2))*d2);
    
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
function plotctrol(Obtained_Pareto,True_Pareto)
if M == 2
    plot(Obtained_Pareto(:,1),Obtained_Pareto(:,2),'o','LineWidth',2,...
        'MarkerEdgeColor','r','MarkerSize',2);
    hold on
    plot(True_Pareto(:,1),True_Pareto(:,2),'k'); 
    title('Optimal Solution Pareto Set using MOSMA');
    legend('MOSMA');
    xlabel('F_1');
    ylabel('F_2');
elseif M == 3
    plot3(Obtained_Pareto(:,1),Obtained_Pareto(:,2),Obtained_Pareto(:,3),'o','LineWidth',2,...
        'MarkerEdgeColor','r','MarkerSize',2);
    hold on
    plot3(Obtained_Pareto(:,1),Obtained_Pareto(:,2),Obtained_Pareto(:,3),'.','LineWidth',2,...
        'MarkerEdgeColor','k','MarkerSize',6);
    title('Optimal Solution Pareto Set using MOSMA');
    legend('MOSMA');
    xlabel('F_1');
    ylabel('F_2');
    zlabel('F_3');
end
