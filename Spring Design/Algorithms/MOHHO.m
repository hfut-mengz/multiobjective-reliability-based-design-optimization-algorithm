function rst=MOHHO(funj,filename,nloop,nsol,nvari,nbit,narchive,a,b);

fun=funj;
%% Problem Definition

% CostFunction=@(x) CostFunction(x);      % Cost Function

nVar=nvari;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables


%% MOPSO Parameters

MaxIt=nloop;           % Maximum Number of Iterations

nPop=nsol;            % Population Size

nRep=50;            % Repository Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

nGrid=7;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate

%% Initialization

empty_Rabbit.Location=[];
empty_Rabbit.Cost=[];
empty_Rabbit.Sol=[];
empty_Rabbit.IsDominated=[];
empty_Rabbit.GridIndex=[];
empty_Rabbit.GridSubIndex=[];

Rabbits=repmat(empty_Rabbit,nPop,1);
X = zeros(nPop, nVar);
Rabbit_Location=zeros(VarSize);
Rabbit_Energy=inf;

for i=1:nPop
    
    Rabbits(i).Location = rand(VarSize).*(VarMax-VarMin)+VarMin; 
    X(i,:) = rand(VarSize).*(VarMax-VarMin)+VarMin; 
%     [Rabbits(i).Cost, Rabbits(i).Sol] = CostFunction(Rabbits(i).Location);
    [Rabbits(i).Cost]=feval(fun,(Rabbits(i).Location)');
    
end

% Determine Domination
Rabbits=DetermineDomination(Rabbits);

rep=Rabbits(~[Rabbits.IsDominated]);

Grid=CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end


%% MOPSO Main Loop

for it=1:MaxIt
    E1=2*(1-(it/MaxIt)); % factor to show the decreaing energy of rabbit    
    for i=1:nPop
        
        leader=SelectLeader(rep,beta);
        
        
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(nPop*rand()+1);
            X_rand = Rabbits(rand_Hawk_index);
            if q<0.5
                % perch based on other family members
                Rabbits(i).Location=X_rand.Location-rand()*abs(X_rand.Location-2*rand()*Rabbits(i).Location);
                X(i,:)=X_rand.Location-rand()*abs(X_rand.Location-2*rand()*Rabbits(i).Location);
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                Rabbits(i).Location=(leader.Location-mean(X))-rand()*((VarMax-VarMin)*rand+VarMin);
                X(i,:)=(leader.Location-mean(X))-rand()*((VarMax-VarMin)*rand+VarMin);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                Rabbits(i).Location=(leader.Location)-Escaping_Energy*abs(leader.Location-Rabbits(i).Location);
                X(i,:)=(leader.Location)-Escaping_Energy*abs(leader.Location-X(i,:));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(leader.Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                Rabbits(i).Location=(leader.Location-Rabbits(i).Location)-Escaping_Energy*abs(Jump_strength*Rabbit_Location-Rabbits(i).Location);
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-X(i,:));
%                 [X1.Cost] = CostFunction(X1.Location);
                [X1.Cost]=feval(fun,X1.Location);
                if Dominates(X1,Rabbits(i))
                    Rabbits(i).Location=X1.Location;
                    Rabbits(i).Cost=X1.Cost;
                    Rabbits(i).Cost=X1.Cost;

                elseif Dominates(Rabbits(i),X1)
                    X2.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-X(i,:))+rand(1,nVar).*Levy(nVar);                    
%                     [X2.Cost] = CostFunction(X2.Location);
                    [X2.Cost]=feval(fun,X2.Location);
                    if Dominates(X2,Rabbits(i))
                        Rabbits(i).Location=X2.Location;
                        Rabbits(i).Cost=X2.Cost;
                        Rabbits(i).Cost=X2.Cost;
                    end
                else
                    if rand<0.5
                        Rabbits(i).Location=X1.Location;
                        Rabbits(i).Cost=X1.Cost;
                        Rabbits(i).Cost=X1.Cost;
%                         Rabbits(i).Sol=X1.Sol;
                    end
                end                                
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-mean(X));
                [X1.Cost] = CostFunction(X1.Location);
                if Dominates(X1,Rabbits(i))
                    Rabbits(i).Location=X1.Location;
                    Rabbits(i).Cost=X1.Cost;
                    Rabbits(i).Cost=X1.Cost;

                elseif Dominates(Rabbits(i),X1)
                    X2.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-mean(X))+rand(1,nVar).*Levy(nVar);                    
%                     [X2.Cost] = CostFunction(X2.Location);
                      [X2.Cost]=feval(fun,X2.Location);
                    if Dominates(X2,Rabbits(i))
                        Rabbits(i).Location=X2.Location;
                        Rabbits(i).Cost=X2.Cost;
                        Rabbits(i).Cost=X2.Cost;
                    end
                else
                    if rand<0.5
                        Rabbits(i).Location=X1.Location;
                        Rabbits(i).Cost=X1.Cost;
                        Rabbits(i).Cost=X1.Cost;
%                         Rabbits(i).Sol=X1.Sol;
                    end
                end                               
            end
        end
        
        Rabbits(i).Location = max(Rabbits(i).Location, VarMin);
        Rabbits(i).Location = min(Rabbits(i).Location, VarMax);
        
%         [Rabbits(i).Cost] = CostFunction(Rabbits(i).Location);
          [Rabbits(i).Cost]=feval(fun,Rabbits(i).Location);
        
        % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            NewSol.Location=Mutate(Rabbits(i).Location,pm,VarMin,VarMax);
%             [NewSol.Cost]=CostFunction(NewSol.Location);
            [NewSol.Cost]=feval(fun,NewSol.Location);
            if Dominates(NewSol,Rabbits(i))
                Rabbits(i).Location=NewSol.Location;
                Rabbits(i).Cost=NewSol.Cost;
%                 Rabbits(i).Sol=NewSol.Sol;

            elseif Dominates(Rabbits(i),NewSol)
                % Do Nothing

            else
                if rand<0.5
                    Rabbits(i).Location=NewSol.Location;
                    Rabbits(i).Cost=NewSol.Cost;
%                     Rabbits(i).Sol=NewSol.Sol;
                end
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         Rabbits(~[Rabbits.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
%     figure(1);
%     PlotCosts(Rabbits,rep);
%     pause(0.01);
    
%     Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w=w*wdamp;
    for i=1:nsol
        dx(i,:)=Rabbits(i).Location;
        fpp(i,:)=(Rabbits(i).Cost)';
    end
    rst.ppareto{it}=dx';
    rst.fpareto{it}=fpp';
    rst.gpareto{it}=[];
    rst.timestamp=datetime('now');
end

end
%% Resluts

% solutions = [];
% costs = [];
% for i=1:numel(rep)
% solutions = cat(1, solutions, rep.Location(i));
% costs = cat(1, costs, rep.Cost(i));
% Sols = cat(1, Sols, rep.Sol(i));
% end
function pop=DetermineDomination(pop)

    nPop=numel(pop);
    
    for i=1:nPop
        pop(i).IsDominated=false;
    end
    
    for i=1:nPop-1
        for j=i+1:nPop
            
            if Dominates(pop(i),pop(j))
               pop(j).IsDominated=true;
            end
            
            if Dominates(pop(j),pop(i))
               pop(i).IsDominated=true;
            end
            
        end
    end

end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function Grid=CreateGrid(pop,nGrid,alpha)

    c=[pop.Cost];
    
    cmin=min(c,[],2);
    cmax=max(c,[],2);
    
    dc=cmax-cmin;
    cmin=cmin-alpha*dc;
    cmax=cmax+alpha*dc;
    
    nObj=size(c,1);
    
    empty_grid.LB=[];
    empty_grid.UB=[];
    Grid=repmat(empty_grid,nObj,1);
    
    for j=1:nObj
        
        cj=linspace(cmin(j),cmax(j),nGrid+1);
        
        Grid(j).LB=[-inf cj];
        Grid(j).UB=[cj +inf];
        
    end

end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function particle=FindGridIndex(particle,Grid)

    nObj=numel(particle.Cost);
    
    nGrid=numel(Grid(1).LB);
    
    particle.GridSubIndex=zeros(1,nObj);
    
    for j=1:nObj
        
        particle.GridSubIndex(j)=...
            find(particle.Cost(j)<Grid(j).UB,1,'first');
        
    end

    particle.GridIndex=particle.GridSubIndex(1);
    for j=2:nObj
        particle.GridIndex=particle.GridIndex-1;
        particle.GridIndex=nGrid*particle.GridIndex;
        particle.GridIndex=particle.GridIndex+particle.GridSubIndex(j);
    end
    
end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function leader=SelectLeader(rep,beta)

    % Grid Index of All Repository Members
    GI=[rep.GridIndex];
    
    % Occupied Cells
    OC=unique(GI);
    
    % Number of Particles in Occupied Cells
    N=zeros(size(OC));
    for k=1:numel(OC)
        N(k)=numel(find(GI==OC(k)));
    end
    
    % Selection Probabilities
    P=exp(-beta*N);
    P=P/sum(P);
    
    % Selected Cell Index
    sci=RouletteWheelSelection(P);
    
    % Selected Cell
    sc=OC(sci);
    
    % Selected Cell Members
    SCM=find(GI==sc);
    
    % Selected Member Index
    smi=randi([1 numel(SCM)]);
    
    % Selected Member
    sm=SCM(smi);
    
    % Leader
    leader=rep(sm);

end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function xnew=Mutate(x,pm,VarMin,VarMax)

    nVar=numel(x);
    j=randi([1 nVar]);

    dx=pm*(VarMax-VarMin);
    
    lb=x(j)-dx;
    if lb<VarMin
        lb=VarMin;
    end
    
    ub=x(j)+dx;
    if ub>VarMax
        ub=VarMax;
    end
    
    xnew=x;
    xnew(j)=unifrnd(lb,ub);

end


function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function b=Dominates(x,y)

    if isstruct(x)
        x=x.Cost;
    end
    
    if isstruct(y)
        y=y.Cost;
    end

    b=all(x<=y) && any(x<y);

end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function rep=DeleteOneRepMemebr(rep,gamma)

    % Grid Index of All Repository Members
    GI=[rep.GridIndex];
    
    % Occupied Cells
    OC=unique(GI);
    
    % Number of Particles in Occupied Cells
    N=zeros(size(OC));
    for k=1:numel(OC)
        N(k)=numel(find(GI==OC(k)));
    end
    
    % Selection Probabilities
    P=exp(gamma*N);
    P=P/sum(P);
    
    % Selected Cell Index
    sci=RouletteWheelSelection(P);
    
    % Selected Cell
    sc=OC(sci);
    
    % Selected Cell Members
    SCM=find(GI==sc);
    
    % Selected Member Index
    smi=randi([1 numel(SCM)]);
    
    % Selected Member
    sm=SCM(smi);
    
    % Delete Selected Member
    rep(sm)=[];

end

%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function i=RouletteWheelSelection(P)

    r=rand;
    
    C=cumsum(P);
    
    i=find(r<=C,1,'first');

end