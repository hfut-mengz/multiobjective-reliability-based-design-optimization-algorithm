% Project Title: A multiobjective mayfly optimization algorithm (MOMA) in MATLAB
%
% Developers: K. Zervoudakis & S. Tsafarakis
%
% Contact Info: kzervoudakis@isc.tuc.gr
%               School of Production Engineering and Management,
%               Technical University of Crete, Chania, Greece
%
% Researchers are allowed to use this code in their research projects.
%
% Please cite as:
% Zervoudakis, K., & Tsafarakis, S. (2020). A mayfly optimization algorithm.
% Computers & Industrial Engineering, 145, 106559.
% https://doi.org/10.1016/j.cie.2020.106559
%%
function rst = MOMA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b) 

%% Problem Definition
% Objective Functions
% ANSWER=listdlg('PromptString','Choose Objective Function','SelectionMode','single', 'ListString', {'1. ZDT', '2. ZDT2', '3. ZDT3'});
% if eq(ANSWER,1); ObjectiveFunction=@(x) ZDT(x); funcname='ZDT';
% elseif eq(ANSWER,2); ObjectiveFunction=@(x) ZDT2(x); funcname='ZDT2';
% elseif eq(ANSWER,3); ObjectiveFunction=@(x) ZDT3(x); funcname='ZDT3';
% else; disp('Terminated'); return
% end
ProblemSize=nvar;       % Decision Variables Size
LowerBound=a(1);             % Decision Variables Lower Bound
UpperBound=b(1);             % Decision Variables Upper Bound
%% Mayfly Parameters
% methname='Mayfly Algorithm';
MaxIt=nloop;         % Maximum Number of Iterations
nPop=nsol; nPopf=nsol;          % Population Size (males and females)
nPareto=nsol;       % Repository Size
g=0.8;                % Inertia Weight
gdamp=1;            % Inertia Weight Damping Ratio
a1=1.0;             % Personal Learning Coefficient
a2=1.5;  a3=1.5;           % Global Learning Coefficient
beta=2;             % Distance sight Coefficient
dance=0.77;          % Mutation Coefficient
dance_damp=0.99;    % Mutation Coefficient Damping Ratio
fl=0.77;                       % Random flight
fl_damp=0.99;
% Mating Parameters
nCrossover=20;  % Number of Parnets (Offsprings)
nMutation=round(0.5*nPop);        % Number of Mutants
mu=0.02;                                % Mutation Rate
% Velocity Limits
VelMax=1*(max(UpperBound)-max(LowerBound))*5; VelMin=-VelMax;
%% Initialization
%run initial
empty_mayfly.Position=[];
empty_mayfly.Velocity=[];
empty_mayfly.Cost=[];
empty_mayfly.Best.Position=[];
empty_mayfly.Best.Cost=[];
empty_mayfly.Rank=[];
empty_mayfly.DominationSet=[];
empty_mayfly.DominatedCount=[];
empty_mayfly.CrowdingDistance=[];
Mayfly=repmat(empty_mayfly,nPop,1);
Mayflyf=repmat(empty_mayfly,nPopf,1);
for i=1:nPop
    % Initialize Male Position
    Mayfly(i).Position=unifrnd(LowerBound,UpperBound,1,ProblemSize);
    
    % Initialize Velocity
    Mayfly(i).Velocity=zeros(1,ProblemSize);
    % Evaluation
%     Mayfly(i).Cost=ObjectiveFunction(Mayfly(i).Position);
    
    [Mayfly(i).Cost,Gi(:,i)]=feval(fun,Mayfly(i).Position');
    Mayfly(i).Cost=fpenal_1(Mayfly(i).Cost,Gi(:,i));
    
    % Update Personal Best
    Mayfly(i).Best.Position=Mayfly(i).Position;
    Mayfly(i).Best.Cost=Mayfly(i).Cost;
    % Initialize female Position
    if i<=nPopf
        Mayflyf(i).Position=unifrnd(LowerBound,UpperBound,1,ProblemSize);
        Mayflyf(i).Velocity=zeros(1,ProblemSize);
%         Mayflyf(i).Cost=ObjectiveFunction(Mayflyf(i).Position);
        
        [Mayflyf(i).Cost,Gi(:,i)]=feval(fun,Mayflyf(i).Position');
        Mayflyf(i).Cost=fpenal_1(Mayflyf(i).Cost,Gi(:,i));
        
        
        Mayflyf(i).Best.Position=Mayflyf(i).Position;
        Mayflyf(i).Best.Cost=Mayflyf(i).Cost;
    end
end
% Merge
Pareto=[Mayfly;Mayflyf];
% Non-Dominated Sorting
[Pareto, F]=ParetoSorting(Pareto);
% Calculate Crowding Distance
Pareto=CalcCD(Pareto,F);
% Sort Population
Pareto=SortSolutions(Pareto);
Pareto=Pareto(F{1});
% Truncate
if numel(Pareto)>nPareto
    Pareto=Pareto(1:nPareto);
end
%% Mayfly Main Loop
for it=1:MaxIt
    for i=1:nPop
        leader=Pareto(randi(size(Pareto,2)));
        % Update Females
        if i<=nPopf
            if Dominates(Mayfly(i),Mayflyf(i))
                rmf=norm(Mayfly(i).Position-Mayflyf(i).Position);
                Mayflyf(i).Velocity = g*Mayflyf(i).Velocity ...
                    +a3*exp(-beta*rmf^2).*(Mayfly(i).Position-Mayflyf(i).Position);
            else
                e=unifrnd(-1,+1,1,ProblemSize);
                Mayflyf(i).Velocity = g*Mayflyf(i).Velocity+fl*(e);
            end
            % Apply Velocity Limits
            Mayflyf(i).Velocity = max(Mayflyf(i).Velocity,VelMin);
            Mayflyf(i).Velocity = min(Mayflyf(i).Velocity,VelMax);
            % Update Position
            Mayflyf(i).Position = Mayflyf(i).Position + Mayflyf(i).Velocity;
            % Velocity Mirror Effect
            IsOutside=(Mayflyf(i).Position<LowerBound | Mayflyf(i).Position>UpperBound);
            Mayflyf(i).Velocity(IsOutside)=-Mayflyf(i).Velocity(IsOutside);
            % Apply Position Limits
            Mayflyf(i).Position = max(Mayflyf(i).Position,LowerBound);
            Mayflyf(i).Position = min(Mayflyf(i).Position,UpperBound);
            % Evaluation
%             Mayflyf(i).Cost=ObjectiveFunction(Mayflyf(i).Position);
            [Mayflyf(i).Cost,Gi(:,i)]=feval(fun,Mayflyf(i).Position');
            Mayflyf(i).Cost=fpenal_1(Mayflyf(i).Cost,Gi(:,i));
            
            
            Mayflyf(i).Best.Position=Mayflyf(i).Position;
            Mayflyf(i).Best.Cost=Mayflyf(i).Cost;
        end
        % Update Males
        % Update Velocity
        if Dominates(leader,Mayfly(i))
            rpbest=norm(Mayfly(i).Best.Position-Mayfly(i).Position);
            rgbest=norm(leader.Position-Mayfly(i).Position);
            Mayfly(i).Velocity = g*Mayfly(i).Velocity ...
                +a1*exp(-beta*rpbest^2).*(Mayfly(i).Best.Position-Mayfly(i).Position) ...
                +a2*exp(-beta*rgbest^2).*(leader.Position-Mayfly(i).Position);
        else
            e=unifrnd(-1,+1,1,ProblemSize);
            Mayfly(i).Velocity = g*Mayfly(i).Velocity+dance*(e);
        end
        % Apply Velocity Limits
        Mayfly(i).Velocity = max(Mayfly(i).Velocity,VelMin);
        Mayfly(i).Velocity = min(Mayfly(i).Velocity,VelMax);
        % Update Position
        Mayfly(i).Position = Mayfly(i).Position + Mayfly(i).Velocity;
        % Velocity Mirror Effect
        IsOutside=(Mayfly(i).Position<LowerBound | Mayfly(i).Position>UpperBound);
        Mayfly(i).Velocity(IsOutside)=-Mayfly(i).Velocity(IsOutside);
        % Apply Position Limits
        Mayfly(i).Position = max(Mayfly(i).Position,LowerBound);
        Mayfly(i).Position = min(Mayfly(i).Position,UpperBound);
        % Evaluation
%         Mayfly(i).Cost=ObjectiveFunction(Mayfly(i).Position);
        [Mayfly(i).Cost,Gi(:,i)]=feval(fun,Mayfly(i).Position');
        Mayfly(i).Cost=fpenal_1(Mayfly(i).Cost,Gi(:,i));
        % Update Personal Best
        if Dominates(Mayfly(i),Mayfly(i).Best)
            Mayfly(i).Best.Position=Mayfly(i).Position;
            Mayfly(i).Best.Cost=Mayfly(i).Cost;
        elseif Dominates(Mayfly(i).Best,Mayfly(i))
            % Do Nothing
        else
            if rand<0.5
                Mayfly(i).Best.Position=Mayfly(i).Position;
                Mayfly(i).Best.Cost=Mayfly(i).Cost;
            end
        end
    end
    % MATE
    popc=repmat(empty_mayfly,nCrossover/2,2);
    for k=1:nCrossover/2
        % Select Parents
        i1=randi(numel(Pareto));
        i2=randi(numel(Pareto));
        %p1=Mayfly(i1).Best;
        %p2=Mayflyf(i2).Best;
        % Apply Crossover
        [popc(k,1).Position, popc(k,2).Position]=Crossover(Pareto(i1).Position,Pareto(i2).Position);
        % Evaluation
        popc(k,1).Position = max(popc(k,1).Position, LowerBound);
        popc(k,1).Position = min(popc(k,1).Position, UpperBound);
%         popc(k,1).Cost=ObjectiveFunction(popc(k,1).Position);
        [ popc(k,1).Cost,Gi(:,i)]=feval(fun,popc(k,1).Position');
         popc(k,1).Cost=fpenal_1( popc(k,1).Cost,Gi(:,i));
        
        % Evaluation
        popc(k,2).Position = max(popc(k,2).Position, LowerBound);
        popc(k,2).Position = min(popc(k,2).Position, UpperBound);
%         popc(k,2).Cost=ObjectiveFunction(popc(k,2).Position);
        [popc(k,2).Cost,Gi(:,i)]=feval(fun,popc(k,2).Position');
         popc(k,2).Cost=fpenal_1( popc(k,2).Cost,Gi(:,i));
        
        popc(k,1).Best.Position = popc(k,1).Position;
        popc(k,1).Best.Cost = popc(k,1).Cost;
        popc(k,1).Velocity= zeros(1,ProblemSize);
        popc(k,2).Best.Position = popc(k,2).Position;
        popc(k,2).Best.Cost = popc(k,2).Cost;
        popc(k,2).Velocity= zeros(1,ProblemSize);
    end
    % break
    popc=popc(:);
    % Mutation
    popm=repmat(empty_mayfly,nMutation,1);
    for k=1:nMutation
        i=randi(numel(Pareto));
        popm(k)=Pareto(i);
        popm(k).Position=Mutate(popm(k).Position,mu,LowerBound,UpperBound);
        % Evaluation
        popm(k).Position = max(popm(k).Position, LowerBound);
        popm(k).Position = min(popm(k).Position, UpperBound);
%         popm(k).Cost=ObjectiveFunction(popm(k).Position);
        
        [popm(k).Cost,Gi(:,i)]=feval(fun,popm(k).Position');
        popm(k).Cost=fpenal_1( popm(k).Cost,Gi(:,i));
    end
    % Create Merged Population
    popc=[popc
        popm]; %#ok
    split=round((nCrossover/2+nMutation)/2);
    males=popc(1:split);
    Mayfly=[Mayfly
        males]; %#ok
    males=popc(split+1:nCrossover/2+nMutation);
    Mayflyf=[Mayflyf
        males]; %#ok
    % SHORT
    % Non-Dominated Sorting
    [Mayfly, F]=ParetoSorting(Mayfly);
    Mayfly=CalcCD(Mayfly,F);
    [Mayfly, F]=SortSolutions(Mayfly);
    [Mayflyf, F]=ParetoSorting(Mayflyf);
    Mayflyf=CalcCD(Mayflyf,F);
    [Mayflyf, F]=SortSolutions(Mayflyf);
    Mayfly=Mayfly(1:nPop);
    Mayflyf=Mayflyf(1:nPopf);
    Pareto=[Pareto
        Mayfly
        Mayflyf]; %#ok
    all=Pareto;
    % Non-Dominated Sorting
    [Pareto, F]=ParetoSorting(Pareto);
    % Calculate Crowding Distance
    Pareto=CalcCD(Pareto,F);
    % Sort Population
    [Pareto, F]=SortSolutions(Pareto);
    % Store F1
    Pareto=Pareto(F{1});
    % Truncate
    if numel(Pareto)>nPareto
        Pareto=Pareto(1:nPareto);
    end
    % Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Number of Solution in repository = ' num2str(numel(Pareto))]);
    % Plot F1 Costs
%     figure(1);
%     PlotCosts(all,Pareto);

% numel(Pareto)
%     str2num(Pareto(1).Position)
%     B=Pareto.Cost;

for j=1:numel(Pareto)
    A.ppareto{j,it}=Pareto(j).Position;
    A.fpareto{j,it}=Pareto(j).Cost;
   [Mayflyf(j).Cost,gi(:,j)]=feval(fun,Pareto(j).Position');
   resutspp(:,j)=Pareto(j).Position';
   resutsfp(:,j)=Pareto(j).Cost';
   resutsgp(:,j)=gi(:,j);
end
% rst.ppareto{it} = mat2cell(A,dim1Dist,...,dimNDist)
    rst.ppareto{it}=resutspp;
    rst.fpareto{it}=resutsfp;
    rst.gpareto{it}=resutsgp;
    rst.timestamp=datetime('now');
    %pause(0.01);
    g=g*gdamp;
    dance = dance*dance_damp;
    fl = fl*fl_damp;
end
%% Results


function [mayflies, F]=SortSolutions(mayflies)
% Sort Based on Crowding Distance
[~, CD]=sort([mayflies.CrowdingDistance],'descend'); mayflies=mayflies(CD);
[~, SO]=sort([mayflies.Rank]); mayflies=mayflies(SO);
Ranking=[mayflies.Rank];
MxRnks=max(Ranking);
F=cell(MxRnks,1);
for r=1:MxRnks
    F{r}=find(Ranking==r);
end


function PlotCosts(mayflies,Pareto)
mayflies_costs=[mayflies.Cost];
plot(mayflies_costs(1,:),mayflies_costs(2,:),'b+','MarkerSize', 4);
hold on;
pareto_costs=[Pareto.Cost];
plot(pareto_costs(1,:),pareto_costs(2,:),'gd','MarkerFaceColor','red','MarkerSize', 10);
grid on;
hold off;


function [mayflies, F]=ParetoSorting(mayflies)
nPop=numel(mayflies);
for i=1:nPop
    mayflies(i).DominationSet=[]; mayflies(i).DominatedCount=0;
end
F{1}=[];
for i=1:nPop
    for j=i+1:nPop
        p=mayflies(i); q=mayflies(j);
        if Dominates(p,q)
            p.DominationSet=[p.DominationSet j];
            q.DominatedCount=q.DominatedCount+1;
        end
        if Dominates(q.Cost,p.Cost)
            q.DominationSet=[q.DominationSet i];
            p.DominatedCount=p.DominatedCount+1;
        end
        mayflies(i)=p; mayflies(j)=q;
    end
    if mayflies(i).DominatedCount==0
        F{1}=[F{1} i]; mayflies(i).Rank=1;
    end
end
k=1;
while true
    Q=[];
    for i=F{k}
        p=mayflies(i);
        for j=p.DominationSet
            q=mayflies(j); q.DominatedCount=q.DominatedCount-1;
            if q.DominatedCount==0
                Q=[Q j]; %#ok
                q.Rank=k+1;
            end
            mayflies(j)=q;
        end
    end
    if isempty(Q)
        break;
    end
    F{k+1}=Q; %#ok
    k=k+1;
end


function y=Mutate(x,mu,LowerBound,UpperBound)
nVar=numel(x); nmu=ceil(mu*nVar);
j=randsample(nVar,nmu);
sigma(1:nVar)=0.1*(UpperBound-LowerBound);
y=x;
y(j)=x(j)+sigma(j)*(randn(size(j))');
y=max(y,LowerBound); y=min(y,UpperBound);


function b=Dominates(x,y)
if isstruct(x); x=x.Cost; end
if isstruct(y); y=y.Cost; end
b=all(x<=y) && any(x<y);


function [off1, off2]=Crossover(x1,x2)
L=unifrnd(0,1,size(x1));
off1=L.*x1+(1-L).*x2;
off2=L.*x2+(1-L).*x1;


function mayflies=CalcCD(mayflies,F)
nF=numel(F);
for k=1:nF
    Objectives=[mayflies(F{k}).Cost];
    Obj=size(Objectives,1); n=numel(F{k});
    d=zeros(n,Obj);
    for j=1:Obj
        [cj, so]=sort(Objectives(j,:)); d(so(1),j)=inf;
        for i=2:n-1
            d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
        end
        d(so(end),j)=inf;
    end
    for i=1:n
        mayflies(F{k}(i)).CrowdingDistance=sum(d(i,:));
    end
end

function [fp] = fpenal_1(f,g)
    % penalty function
    if max(g)>0
        fp=f+100*max(g);
    else
        fp=f;
    end
