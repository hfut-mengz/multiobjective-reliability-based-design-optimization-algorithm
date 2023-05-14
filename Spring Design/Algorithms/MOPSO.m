function rst = MOPSO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
%_________________________________________________________________________________
%  Multi-objective Particle Swarm Optimization (MOPSO) source codes version 1.0
%
%  Developed in MATLAB R2016a
%
%  Author and programmer: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%                 seyedali.mirjalili@griffithuni.edu.au
%
%       Homepage: http://www.alimirjalili.com
%
%____________________________________________________________________________________
rand('state',sum(100*clock));
Vmax=10;
wMax=0.9;
wMin=0.4;
c1=1.5;
c2=1.5;

% Change these details with respect to your problem%%%%%%%%%%%%%%
dim=nvar;
lb=a';
ub=b';
obj_no=numel(feval(fun,a));

if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Max_it=nloop;
N=nsol;
ArchiveMaxSize=narchive;
vel=zeros(dim,N);
Archive_X=zeros(ArchiveMaxSize,dim);
[ff,gg]=feval(fun,Archive_X(1,:)');
Archive_F=ones(ArchiveMaxSize,obj_no)*inf;
Archive_G=ones(ArchiveMaxSize,numel(gg))*inf;

Archive_member_no=0;

r=(ub-lb)/2;
V_max=(ub(1)-lb(1))/10;

Food_fitness=inf*ones(1,obj_no);
Food_position=zeros(dim,1);
Food_G=ones(1,numel(gg))*inf;

Salps_X=initialization(N,dim,ub,lb);
Salps_fitness=zeros(N,numel(ff));

V=initialization(N,dim,ub,lb);
T=0;

position_history=zeros(N,Max_it,dim);

for T=1:Max_it
    
    for i=1:N
        % Objective Function
        [fi,gi]=feval(fun,Salps_X(:,i));
        % Penalty Function
        if max(gi)>0
            Salps_fitness(i,:)=fpenal_1(fi',gi');
        else
            Salps_fitness(i,:)=fi';
        end
        Salps_G(i,:)=gi';
        %
        %
        if dominates(Salps_fitness(i,:),Salps_G(i,:),Food_fitness,Food_G)
            Food_fitness=Salps_fitness(i,:);
            Food_position=Salps_X(:,i);
            Food_G=gi';
        end
    end
    
    [Archive_X, Archive_F, Archive_G, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, Salps_X, Salps_fitness, Salps_G, Archive_member_no);
    
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
        [Archive_X, Archive_F, Archive_G, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
    end
    
    Archive_mem_ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no);
    % Archive_mem_ranks
    % Chose the archive member in the least population area as food`
    % to improve coverage
    index=RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
    Food_fitness=Archive_F(index,:);
    Food_position=Archive_X(index,:)';
    
    % PSO updating mechanism
    w=wMax-T*((wMax-wMin)/Max_it);
    index=0;
    neighbours_no=0;
    for i=1:N
        for j = 1:dim
            if i > size(Archive_X,1)
                RJ = size(Archive_X,1);
            else
                RJ = i;
            end
            vel(j,i) = w*vel(j,i)+c1*rand*(Archive_X(RJ,j)-Salps_X(j,i))+...
                c2*rand*(Food_position(j)-Salps_X(j,i));   % velocity updating
            if vel(j,i) > Vmax
                vel(j,i) = Vmax;
            end
            if vel(j,i) < -Vmax
                vel(j,i) = -Vmax;
            end
            Salps_X(j,i) = Salps_X(j,i)+vel(j,i);       % position updating
        end
        Flag4ub=Salps_X(:,i)>ub';
        Flag4lb=Salps_X(:,i)<lb';
        Salps_X(:,i)=(Salps_X(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
    end
    
    %     display(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
    
    % Save results
    [ppareto,fpareto,gpareto,~]=pbil_selection0([],[],[],Archive_X',Archive_F',Archive_G',narchive);
    rst.ppareto{T}=ppareto;
    rst.fpareto{T}=fpareto;
    rst.gpareto{T}=gpareto;
    rst.timestamp=datetime('now');
end


ppareto=Archive_X';
for i=1:size(Archive_X,1)
    [fpareto(:,i),gpareto(:,i)]=feval(fun,ppareto(:,i));
end

end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function o=dominates(f1,g1,f2,g2)
if numel(f1)>0 && numel(g1)>0 && numel(f2)>0 && numel(g2)>0
    if max(g1)<=0 && max(g2)<=0
        o=all(f1<=f2) && any(f1<f2);
    elseif max(g1)<=0 && max(g2)>0
        o=1;
    elseif max(g1)>0 && max(g2)<=0
        o=0;
    elseif max(g1)>0 && max(g2)>0
        if max(g1)>=max(g2)
            o=0;
        else
            o=1;
        end
    else
        error('unknown error, please check');
    end
elseif numel(f1)>0 && numel(g1)>0 && numel(f2)==0 && numel(g2)==0
    o=1;
elseif numel(f1)==0 && numel(g1)==0 && numel(f2)>0 && numel(g2)>0
else
    o=0;
end
end

function [Archive_X_Chopped, Archive_F_Chopped, Archive_G_Chopped, Archive_mem_ranks_updated, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_G, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize)
for i=1:size(Archive_F,1)-ArchiveMaxSize
    index=RouletteWheelSelection(Archive_mem_ranks);
    %[value index]=min(Archive_mem_ranks);
    
    if index==-1
        index=ceil(size(Archive_F,1)*rand(size(Archive_F,1)-ArchiveMaxSize,1));
    end
    
    Archive_X=[Archive_X(1:index-1,:) ; Archive_X(index+1:Archive_member_no,:)];
    Archive_F=[Archive_F(1:index-1,:) ; Archive_F(index+1:Archive_member_no,:)];
    Archive_G=[Archive_G(1:index-1,:) ; Archive_G(index+1:Archive_member_no,:)];
    Archive_mem_ranks=[Archive_mem_ranks(1:index-1) Archive_mem_ranks(index+1:Archive_member_no)];
    Archive_member_no=Archive_member_no-1;
end

Archive_X_Chopped=Archive_X;
Archive_F_Chopped=Archive_F;
Archive_G_Chopped=Archive_G;
Archive_mem_ranks_updated=Archive_mem_ranks;
end

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    ub_new=ones(1,dim)*ub;
    lb_new=ones(1,dim)*lb;
else
    ub_new=ub;
    lb_new=lb;
end

% If each variable has a different lb and ub
for i=1:dim
    ub_i=ub_new(i);
    lb_i=lb_new(i);
    Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
end

Positions=Positions';
end

function ranks=RankingProcess(Archive_F, Archive_G, ArchiveMaxSize, obj_no)

my_min=min(Archive_F);
my_max=max(Archive_F);

if size(Archive_F,1)==1
    my_min=Archive_F;
    my_max=Archive_F;
end

r=(my_max-my_min)/(20);

ranks=zeros(1,size(Archive_F,1));

for i=1:size(Archive_F,1)
    ranks(i)=0;
    for j=1:size(Archive_F,1)
        flag=0; % a flag to see if the point is in the neoghbourhood in all dimensions.
        for k=1:obj_no
            %
            %
            F1=fpenal_1(Archive_F(j,k),Archive_G(j,:));
            F2=fpenal_1(Archive_F(i,k),Archive_G(i,:));
            if (abs(F1-F2)<r(k))
                flag=flag+1;
            end
            %
            %
        end
        if flag==obj_no
            ranks(i)=ranks(i)+1;
        end
    end
end
end

function choice = RouletteWheelSelection(weights)
accumulation = cumsum(weights);
p = rand() * accumulation(end);
chosen_index = -1;
for index = 1 : length(accumulation)
    if (accumulation(index) > p)
        chosen_index = index;
        break;
    end
end
choice = chosen_index;
end

function [Archive_X_updated, Archive_F_updated, Archive_G_updated, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Archive_G, Particles_X, Particles_F, Particles_G, Archive_member_no)
[Archive_X_temp,uind]=unique([Archive_X ; Particles_X'],'rows');
Archive_F_temp=[Archive_F ; Particles_F];
Archive_F_temp=Archive_F_temp(uind,:);
Archive_G_temp=[Archive_G ; Particles_G];
Archive_G_temp=Archive_G_temp(uind,:);

o=zeros(1,size(Archive_F_temp,1));

for i=1:size(Archive_F_temp,1)
    o(i)=0;
    for j=1:i-1
        if any(Archive_F_temp(i,:) ~= Archive_F_temp(j,:))
            if dominates(Archive_F_temp(i,:),Archive_G_temp(i,:),Archive_F_temp(j,:),Archive_G_temp(j,:))
                o(j)=1;
            elseif dominates(Archive_F_temp(j,:),Archive_G_temp(j,:),Archive_F_temp(i,:),Archive_G_temp(i,:))
                o(i)=1;
                break;
            end
        else
            o(j)=0;
            o(i)=0;
        end
    end
end

Archive_member_no=0;
index=0;
for i=1:size(Archive_X_temp,1)
    if o(i)==0
        Archive_member_no=Archive_member_no+1;
        Archive_X_updated(Archive_member_no,:)=Archive_X_temp(i,:);
        Archive_F_updated(Archive_member_no,:)=Archive_F_temp(i,:);
        Archive_G_updated(Archive_member_no,:)=Archive_G_temp(i,:);
    else
        index=index+1;
        %         dominated_X(index,:)=Archive_X_temp(i,:);
        %         dominated_F(index,:)=Archive_F_temp(i,:);
    end
end
end


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
end

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
end

function [fp] = fpenal_1(f,g)
% penalty function
if max(g)>0
    fp=f+100*max(g);
else
    fp=f;
end
end