
function rst = NSMFO(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
tic
fobj=fun;
LB=a';
UB=b';
D=nvar;
Max_iteration=nloop;
SearchAgents_no=nsol;
M=numel(feval(fun,LB'));
%% Non Sorted Moth-flame Optimization Algorithm (NSMFO)
% NSMFO is developed by Pradeep Jangir
% f - optimal fitness
% X - optimal solution
% D  Dimensional of problem at hand   
% M Number of objective function
% Moth_pos is a matrix consists of all individuals
% SearchAgents_no is number of individual in Moth_possystem
% LB lower boundary constraint
% UB upper boundary constraint
%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.
chromosome = initialize_variables(SearchAgents_no, M, D, LB, UB,fun);
%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
intermediate_chromosome = non_domination_sort_mod(chromosome, M, D);
%% Perform Selection
% Once the intermediate population is sorted only the best solution is
% selected based on it rank and crowding distance. Each front is filled in
% ascending order until the addition of population size is reached. The
% last front is included in the population based on the individuals with
% least crowding distance
% Select NP fittest solutions using non dominated and crowding distance
% sorting and store in population
Pop = replace_chromosome(intermediate_chromosome, M,D,SearchAgents_no);
%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.
%% Algorithm Variables
K = D+M;
Moth_pos = Pop(:,1:K+1);
Moth_pos_ad = zeros(SearchAgents_no,K);
%% Optimization Circle
Iteration = 1;
while Iteration<=Max_iteration % for each generation    
    for i = 1:SearchAgents_no   %  (Moth for each individual)
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        SF=round(1+rand); %% Scaling factor to perform best coverage in MFO
        % randomly select the best organism from first non-dominated front of Moth_pos
        Ffront = Moth_pos((find(Moth_pos(:,K+1)==1)),:); % first front
        ri =  floor(size(Ffront,1)*rand())+1; % ri is random index
        sorted_population = Moth_pos(ri,1:D);
     
         Flame_no=round(SearchAgents_no-Iteration*((SearchAgents_no-1)/Max_iteration));
         sorted_population1 = Moth_pos(Flame_no,1:D);
         
        % Calculate new solution 
        Moth_posNew1 = Moth_pos(i,1:D)+rand(1,D).*(sorted_population-SF.*Moth_pos(i,1:D));
        % Handling constraints
        Moth_posNew1 = bound(Moth_posNew1(:,1:D),UB,LB); 
        % Evaluate function at Moth_posNew1 
%         Moth_posNew1(:,D + 1: K) = ZDT1(Moth_posNew1(:,1:D));
        [Moth_posNew1(:,D + 1: K),gi]=feval(fobj,Moth_posNew1(:,1:D)');
%         [Moth_posNew1(:,D + 1: K),gi]=feval(fobj,Moth_posNew1(:,1:D)');
        Moth_posNew1(:,D + 1: K)=fpenal_1(Moth_posNew1(:,D + 1: K),gi);
        % For the first trail Moth_posNew1
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Moth_posNew1(:,D+k)<Moth_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Moth_posNew1(:,D+k)== Moth_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Moth_posNew1 dominates
            % target vector Xi. Replace Xi by Moth_posNew1 in current Moth_possystem and
            % add Xi to advanced population 2
            Moth_pos_ad(i,1:K) = Moth_pos(i,1:K); % add Xi to advanced Moth_pos
            Moth_pos(i,1:K) = Moth_posNew1(:,1:K); % replace Xi by Moth_posNew1            
        else % else Add Xi (trial vector) to advanced Moth_pos
            Moth_pos_ad(i,1:K)= Moth_posNew1;
        end % end if
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
                dom_more = dom_more +1;           
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Moth_posNew1 dominates
            Moth_pos_ad(j,1:K) = Moth_pos(j,1:K); % add Xi to advanced Moth_pos
        end % end if
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end      
        for i=1:size(Moth_pos,1)
            if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
        % Calculate new solution 
        a=-1+Iteration*((-1)/Max_iteration);
        b=1;
        t=(a-1)*rand+1;                
        Moth_posNew1 = Moth_pos(i,1:D)+abs(sorted_population-Moth_pos(j,1:D))*exp(b.*t).*cos(t.*2*pi)+sorted_population;
            end
            if i>Flame_no % Update the position of the moth with respect to its corresponsing flame
        a=-1+Iteration*((-1)/Max_iteration);
        b=1;
        t=(a-1)*rand+1;                
        Moth_posNew1 = Moth_pos(i,1:D)+abs(sorted_population-Moth_pos(j,1:D))*exp(b.*t).*cos(t.*2*pi)+sorted_population1;
            end
        end  
        Moth_posNew1 = bound(Moth_posNew1(:,1:D),UB,LB);
%         Moth_posNew2=fpenal_1Moth_posNew1
%         Moth_posNew1(:,D + 1: K) = evaluate_objective(Moth_posNew1(:,1:D));
        [Moth_posNew1(:,D + 1: K),gi]=feval(fobj,Moth_posNew1(:,1:D)');
        Moth_posNew1(:,D + 1: K)=fpenal_1(Moth_posNew1(:,D + 1: K),gi);
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Moth_posNew1(:,D+k)<Moth_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Moth_posNew1(:,D+k)== Moth_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Moth_posNew1 dominates
            % target vector Xi. Replace Xi by Moth_posNew1 in current Moth_possystem and
            % add Xi to advanced population
            Moth_pos_ad(i,1:K) = Moth_pos(i,1:K); % add Xi to advanced Moth_pos
            Moth_pos(i,1:K) = Moth_posNew1(:,1:K); % replace Xi by Moth_posNew1            
        else % else Add Xi (trial vector) to advanced Moth_pos
            Moth_pos_ad(i,1:K)= Moth_posNew1;
        end % end if
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        parasiteVector=Moth_pos(i,1:D);
        seed=randperm(D);
        pick=seed(1:ceil(rand*D));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(UB(pick)-LB(pick))+LB(pick);        
        % Evaluate the Parasite Vector
%         parasiteVector(:,D + 1: K) = evaluate_objective(parasiteVector(:,1:D));
        [parasiteVector(:,D + 1: K),gi]=feval(fobj,Moth_posNew1(:,1:D)');
        parasiteVector(:,D + 1: K)=fpenal_1(parasiteVector(:,D + 1: K),gi);
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (parasiteVector(:,D+k)<Moth_pos(j,D+k))
                dom_less = dom_less + 1;
            elseif (parasiteVector(:,D+k)== Moth_pos(j,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Moth_posNew1 dominates
            % target vector Xi. Replace Xi by Moth_posNew1 in current Moth_possystem and
            % add Xi to advanced population
            Moth_pos_ad(j,1:K) = Moth_pos(j,1:K); % add Xi to advanced Moth_pos
            Moth_pos(j,1:K) = parasiteVector(:,1:K); % replace Xi by Moth_posNew1            
        else % else Add Xi (trial vector) to advanced Moth_pos
            Moth_pos_ad(j,1:K)= parasiteVector;
        end % end if 
    end % end for i
%     if rem(Iteration, ishow) == 0
%         fprintf('Generation: %d\n', Iteration);        
%     end
    Moth_pos_com = [Moth_pos(:,1:K) ; Moth_pos_ad];
    intermediate_Moth_pos = non_domination_sort_mod(Moth_pos_com, M, D);
    Pop  = replace_chromosome(intermediate_Moth_pos, M,D,SearchAgents_no);
    Moth_pos=Pop(:,1:K+1); %
    rst.ppareto{Iteration}=Moth_pos(:,1:D)'; %设计变量
    rst.fpareto{Iteration}=Moth_pos(:,D+1:K)'; %目标函数
    rst.gpareto{Iteration}=gi; %约束条件,目前要修改
    rst.timestamp=datetime('now');
    Iteration = Iteration+1;
end 
f= Moth_pos;



%% Cited from NSGA-II All rights reserved.
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
    [f(i,D + 1: K),gi] = feval(fun,f(i,1:D)');
    f(i,D + 1: K)=fpenal_1(f(i,D + 1: K),gi);
end % end for i



%% Cited from NSGA-II All rights reserved.
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
%% References
% [1] *Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan*, |A Fast
% Elitist Multiobjective Genetic Algorithm: NSGA-II|, IEEE Transactions on 
% Evolutionary Computation 6 (2002), no. 2, 182 ~ 197.
%
% [2] *N. Srinivas and Kalyanmoy Deb*, |Multiobjective Optimization Using 
% Nondominated Sorting in Genetic Algorithms|, Evolutionary Computation 2 
% (1994), no. 3, 221 ~ 248.



%% Cited from NSGA-II All rights reserved.
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



% Check the boundary limit
function a=bound(a,ub,lb)
a(a>ub)=ub(a>ub); a(a<lb)=lb(a<lb);



function [fp] = fpenal_1(f,g)
    % penalty function
    if max(g)>0
        fp=f+100*max(g);
    else
        fp=f;
    end

    

