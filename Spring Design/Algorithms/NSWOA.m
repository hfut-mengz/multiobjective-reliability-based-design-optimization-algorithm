function rst = NSWOA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)

%% Non Sorted Whale Optimization Algorithm (NSWOA)
% NSWOA is developed by Pradeep Jangir
% f - optimal fitness
% X - optimal solution
% D  Dimensional of problem at hand   
% M Number of objective function
% Whale_pos is a matrix consists of all individuals
% SearchAgents_no is number of individual in Whale_possystem
% LB lower boundary constraint
% UB upper boundary constraint
%% Algorithm Variables
D=nvar;
u0=(a+b)/2;
[f0,g0]=feval(fun,u0);
M=length(f0);
LB=a;
UB=b;
SearchAgents_no=nsol;
Max_iteration=nloop;

K = D+M;

chromosome = initialize_variables(SearchAgents_no, M, D, LB, UB,fun);
intermediate_chromosome = non_domination_sort_mod(chromosome, M, D);
Population = replace_chromosome(intermediate_chromosome, M,D,SearchAgents_no);
Pop=Population;


Whale_pos = Pop(:,1:K+1);
Whale_pos_ad = zeros(SearchAgents_no,K);
%% Optimization Circle
Iteration = 1;
while Iteration<=Max_iteration % for each generation    
    for i = 1:SearchAgents_no   %  (Moth for each individual)
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        SF=round(1+rand); %% Scaling factor to perform best coverage in MFO
        % randomly select the best organism from first non-dominated front of Whale_pos
        Ffront = Whale_pos((find(Whale_pos(:,K+1)==1)),:); % first front
        ri =  floor(size(Ffront,1)*rand())+1; % ri is random index
        sorted_population = Whale_pos(ri,1:D);
        % Calculate new solution 
        Whale_posNew1 = Whale_pos(i,1:D)+rand(1,D).*(sorted_population-SF.*Whale_pos(i,1:D));
        % Handling constraints
        Whale_posNew1 = bound(Whale_posNew1(:,1:D),UB',LB'); 
        % Evaluate function at Whale_posNew1 
%         Whale_posNew1(:,D + 1: K) = evaluate_objective(Whale_posNew1(:,1:D));
        
           [ Whale_posNew1(:,D + 1: K),G(:,i)] = feval(fun,Whale_posNew1(:,1:D)');
           Whale_posNew1(:,D + 1: K)=fpenal_1(Whale_posNew1(:,D + 1: K),G(:,i));
        % For the first trail Whale_posNew1
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Whale_posNew1(:,D+k)<Whale_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Whale_posNew1(:,D+k)== Whale_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            % target vector Xi. Replace Xi by Whale_posNew1 in current Whale_possystem and
            % add Xi to advanced population 2
            Whale_pos_ad(i,1:K) = Whale_pos(i,1:K); % add Xi to advanced Whale_pos
            Whale_pos(i,1:K) = Whale_posNew1(:,1:K); % replace Xi by Whale_posNew1            
        else % else Add Xi (trial vector) to advanced Whale_pos
            Whale_pos_ad(i,1:K)= Whale_posNew1;
        end % end if
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
                dom_more = dom_more +1;           
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            Whale_pos_ad(j,1:K) = Whale_pos(j,1:K); % add Xi to advanced Whale_pos
        end % end if
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end 
        a=2-Iteration*((2)/Max_iteration ); % a decreases linearly fron 2 to 0 in Eq. (2.3)
        a2=-1+Iteration*((-1)/Max_iteration );
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]       
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper 
        b=1;               %  parameters in Eq. (2.5)
        t=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        p = rand();        % p in Eq. (2.6) 
            if p<0.5 % Update the position of the moth with respect to its corresponsing flame
        % Calculate new solution       
        X_rand = sorted_population;
        Whale_posNew1 = Whale_pos(i,1:D)+X_rand-A.*abs(C*X_rand-Whale_pos(j,1:D));
            elseif p>=0.5  
        Whale_posNew1 = Whale_pos(i,1:D)+abs(sorted_population-Whale_pos(j,1:D))*exp(b.*t).*cos(t.*2*pi)+sorted_population;
            end
        Whale_posNew1 = bound(Whale_posNew1(:,1:D),UB',LB');
%         Whale_posNew1(:,D + 1: K) = evaluate_objective(Whale_posNew1(:,1:D));
        
        [Whale_posNew1(:,D + 1: K),Gi(:,i)]=feval(fun,Whale_posNew1(:,1:D)');
        Whale_posNew1(:,D + 1: K)=fpenal_1(Whale_posNew1(:,D + 1: K),Gi(:,i));
        
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Whale_posNew1(:,D+k)<Whale_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Whale_posNew1(:,D+k)== Whale_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            % target vector Xi. Replace Xi by Whale_posNew1 in current Whale_possystem and
            % add Xi to advanced population
            Whale_pos_ad(i,1:K) = Whale_pos(i,1:K); % add Xi to advanced Whale_pos
            Whale_pos(i,1:K) = Whale_posNew1(:,1:K); % replace Xi by Whale_posNew1            
        else % else Add Xi (trial vector) to advanced Whale_pos
            Whale_pos_ad(i,1:K)= Whale_posNew1;
        end % end if
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        parasiteVector=Whale_pos(i,1:D);
        seed=randperm(D);
        pick=seed(1:ceil(rand*D));  % select random dimension
%         aa=rand(1,length(pick));
%         bb=(UB(pick)-LB(pick))';
%         cc=LB(pick)';
        parasiteVector(:,pick)=rand(1,length(pick)).*(UB(pick)'-LB(pick)')+LB(pick)';        
        % Evaluate the Parasite Vector
%         parasiteVector(:,D + 1: K) = evaluate_objective(parasiteVector(:,1:D));
        [parasiteVector(:,D + 1: K),Gi(:,i)]=feval(fun,parasiteVector(:,1:D)');
        parasiteVector(:,D + 1: K)=fpenal_1(parasiteVector(:,D + 1: K),Gi(:,i));
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (parasiteVector(:,D+k)<Whale_pos(j,D+k))
                dom_less = dom_less + 1;
            elseif (parasiteVector(:,D+k)== Whale_pos(j,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            % target vector Xi. Replace Xi by Whale_posNew1 in current Whale_possystem and
            % add Xi to advanced population
            Whale_pos_ad(j,1:K) = Whale_pos(j,1:K); % add Xi to advanced Whale_pos
            Whale_pos(j,1:K) = parasiteVector(:,1:K); % replace Xi by Whale_posNew1            
        else % else Add Xi (trial vector) to advanced Whale_pos
            Whale_pos_ad(j,1:K)= parasiteVector;
        end % end if 
    end % end for i
    if rem(Iteration, 10) == 0
        fprintf('Generation: %d\n', Iteration);        
    end
    Whale_pos_com = [Whale_pos(:,1:K) ; Whale_pos_ad];
    intermediate_Whale_pos = non_domination_sort_mod(Whale_pos_com, M, D);
    Pop  = replace_chromosome(intermediate_Whale_pos, M,D,SearchAgents_no);
    Whale_pos=Pop(:,1:K+1); %    
    
    
    for o=1:size(Pop,1)
        [P(o,:),Gii(:,o)]=feval(fun,Pop(o,1:D)');
    end
    
    rst.ppareto{Iteration }=Pop(:,1:D)';
    rst.fpareto{Iteration }=P';
    rst.gpareto{Iteration }=Gii;
    rst.timestamp=datetime('now');
    
    
    Iteration = Iteration+1;
end 


% Check the boundary limit
function a=bound(a,ub,lb)
a(a>ub)=ub(a>ub); a(a<lb)=lb(a<lb);


function f  = replace_chromosome(intermediate_chromosome, M,D,NP)

%% function f  = replace_chromosome(intermediate_chromosome,M,D,NP)
% This function replaces the chromosomes based on rank and crowding
% distance. Initially until the population size is reached each front is
% added one by one until addition of a complete front which results in
% exceeding the population size. At this point the chromosomes in that
% front is added subsequently to the population based on crowding distance.

[~,m]=size(intermediate_chromosome);
f=zeros(NP,m);

% Now sort the individuals based on the index
sorted_chromosome = sortrows(intermediate_chromosome,M + D + 1);

% Find the maximum rank in the current population
max_rank = max(intermediate_chromosome(:,M + D + 1));

% Start adding each front based on rank and crowing distance until the
% whole population is filled.
previous_index = 0;
for i = 1 : max_rank
    % Get the index for current rank i.e the last the last element in the
    % sorted_chromosome with rank i. 
    current_index = find(sorted_chromosome(:,M + D + 1) == i, 1, 'last' );
    % Check to see if the population is filled if all the individuals with
    % rank i is added to the population. 
    if current_index > NP
        % If so then find the number of individuals with in with current
        % rank i.
        remaining = NP - previous_index;
        % Get information about the individuals in the current rank i.
        temp_pop = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        % Sort the individuals with rank i in the descending order based on
        % the crowding distance.
        [~,temp_sort_index] = ...
            sort(temp_pop(:, M + D + 2),'descend');
        % Start filling individuals into the population in descending order
        % until the population is filled.
        for j = 1 : remaining
            f(previous_index + j,:) = temp_pop(temp_sort_index(j),:);
        end
        return;
    elseif current_index < NP
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
    else
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        return;
    end % end if current_index
    % Get the index for the last added individual.
    previous_index = current_index;
end


function f = non_domination_sort_mod(x, M, D)

%% function f = non_domination_sort_mod(x, M, D)
% This function sort the current popultion based on non-domination. All the
% individuals in the first front are given a rank of 1, the second front
% individuals are assigned rank 2 and so on. After assigning the rank the
% crowding in each front is calculated.

[N, ~] = size(x);

% Initialize the front number to 1.
front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
F(front).f = [];
individual = [];

%% Non-Dominated sort. 
% The initialized population is sorted based on non-domination. The fast
% sort algorithm [1] is described as below for each

% ? for each individual p in main population P do the following
%   ? Initialize Sp = []. This set would contain all the individuals that is
%     being dominated by p.
%   ? Initialize np = 0. This would be the number of individuals that domi-
%     nate p.
%   ? for each individual q in P
%       * if p dominated q then
%           ? add q to the set Sp i.e. Sp = Sp ? {q}
%       * else if q dominates p then
%           ? increment the domination counter for p i.e. np = np + 1
%   ? if np = 0 i.e. no individuals dominate p then p belongs to the first
%     front; Set rank of individual p to one i.e prank = 1. Update the first
%     front set by adding p to front one i.e F1 = F1 + {p}
% ? This is carried out for all the individuals in main population P.
% ? Initialize the front counter to one. i = 1
% ? following is carried out while the ith front is nonempty i.e. Fi != []
%   ? Q = []. The set for storing the individuals for (i + 1)th front.
%   ? for each individual p in front Fi
%       * for each individual q in Sp (Sp is the set of individuals
%         dominated by p)
%           ? nq = nq-1, decrement the domination count for individual q.
%           ? if nq = 0 then none of the individuals in the subsequent
%             fronts would dominate q. Hence set qrank = i + 1. Update
%             the set Q with individual q i.e. Q = Q + q.
%   ? Increment the front counter by one.
%   ? Now the set Q is the next front and hence Fi = Q.
%
% This algorithm is better than the original NSGA ([2]) since it utilize
% the informatoion about the set that an individual dominate (Sp) and
% number of individuals that dominate the individual (np).

%
for i = 1 : N
    % Number of individuals that dominate this individual
    individual(i).n = 0;
    % Individuals which this individual dominate
    individual(i).p = [];
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,D + k) < x(j,D + k))
                dom_less = dom_less + 1;
            elseif (x(i,D + k) == x(j,D + k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).p = [individual(i).p j];
        end
    end   % end for j (N)
    if individual(i).n == 0
        x(i,M + D + 1) = 1;
        F(front).f = [F(front).f i];
    end
end % end for i (N)
% Find the subsequent fronts
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)
       if ~isempty(individual(F(front).f(i)).p)
        	for j = 1 : length(individual(F(front).f(i)).p)
            	individual(individual(F(front).f(i)).p(j)).n = ...
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0
               		x(individual(F(front).f(i)).p(j),M + D + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

sorted_based_on_front = sortrows(x,M+D+1); % sort follow front
current_index = 0;

%% Crowding distance
%The crowing distance is calculated as below
% ? For each front Fi, n is the number of individuals.
%   ? initialize the distance to be zero for all the individuals i.e. Fi(dj ) = 0,
%     where j corresponds to the jth individual in front Fi.
%   ? for each objective function m
%       * Sort the individuals in front Fi based on objective m i.e. I =
%         sort(Fi,m).
%       * Assign infinite distance to boundary values for each individual
%         in Fi i.e. I(d1) = Inf and I(dn) = Inf
%       * for k = 2 to (n - 1)
%           ? I(dk) = I(dk) + (I(k + 1).m - I(k - 1).m)/fmax(m) - fmin(m)
%           ? I(k).m is the value of the mth objective function of the kth
%             individual in I

% Find the crowding distance for each individual in each front
for front = 1 : (length(F) - 1)
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);
    end
    current_index = current_index + i;            
    for i = 1 : M
        [sorted_based_on_objective, index_of_objectives] = sortrows(y,D + i);       
        f_max = ...
            sorted_based_on_objective(length(index_of_objectives), D + i);
        f_min = sorted_based_on_objective(1, D + i);
        y(index_of_objectives(length(index_of_objectives)),M + D + 1 + i)...
            = Inf;
        y(index_of_objectives(1),M + D + 1 + i) = Inf;
         for j = 2 : length(index_of_objectives) - 1
            next_obj  = sorted_based_on_objective(j + 1,D + i);
            previous_obj  = sorted_based_on_objective(j - 1,D + i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + D + 1 + i) = Inf;
            else
                y(index_of_objectives(j),M + D + 1 + i) = ...
                     (next_obj - previous_obj)/(f_max - f_min);
            end
         end % end for j
    end % end for i
    distance = [];
    distance(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        distance(:,1) = distance(:,1) + y(:,M + D + 1 + i);
    end
    y(:,M + D + 2) = distance;
    y = y(:,1 : M + D + 2);
    z(previous_index:current_index,:) = y;
end
f = z();


function f = initialize_variables(NP, M, D, LB, UB,fun)

%% function f = initialize_variables(N, M, D, LB, UB) 
% This function initializes the population. Each individual has the
% following at this stage
%       * set of decision variables
%       * objective function values
% 
% where,
% NP - Population size
% M - Number of objective functions
% D - Number of decision variables
% min_range - A vector of decimal values which indicate the minimum value
% for each decision variable.
% max_range - Vector of maximum possible values for decision variables.

min = LB;
max = UB;

% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.


K = M + D;
f=zeros(NP,K);

%% Initialize each individual in population
% For each chromosome perform the following (N is the population size)
for i = 1 : NP
    % Initialize the decision variables based on the minimum and maximum
    % possible values. V is the number of decision variable. A random
    % number is picked between the minimum and maximum possible values for
    % the each decision variable.
    for j = 1 : D
        f(i,j) = min(j) + (max(j) - min(j))*rand(1);
    end % end for j
    % For ease of computation and handling data the chromosome also has the
    % vlaue of the objective function concatenated at the end. The elements
    % D + 1 to K has the objective function valued. 
    % The function evaluate_objective takes one individual at a time,
    % infact only the decision variables are passed to the function along
    % with information about the number of objective functions which are
    % processed and returns the value for the objective functions. These
    % values are now stored at the end of the individual itself.
%     f(i,D + 1: K) = evaluate_objective(f(i,1:D));
    [f(i,D + 1: K),Gi(:,i)]=feval(fun,f(i,1:D)');
    f(i,D + 1: K)=fpenal_1(f(i,D + 1: K),Gi(:,i));
end % end for i



function [fp] = fpenal_1(f,g)
    % penalty function
    if max(g)>0
        fp=f+100*max(g);
    else
        fp=f;
    end
