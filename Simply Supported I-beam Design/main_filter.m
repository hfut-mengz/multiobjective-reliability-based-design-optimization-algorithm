clc;clear;
filename={'rst_001_005.mat';'rst_001_006.mat';'rst_001_015.mat'};% ;'rst_001_013.mat'
% filename={'rst_001_015.mat'};
%%%%%%%%%%%%计算约束所在的点
for i=1:numel(filename)
    fileiter=filename{i};
    load(fileiter);
    [nrst]=length(rst);
    
    for j=1:nrst
        xtotal=rst(j).ppareto;
        [maxiter]=length(xtotal);
        for k=1:maxiter
            xiter=xtotal{k};
            [nd,np]=size(xiter);
            for t=1:np
                x=xiter(:,t);
                [f,g] = ftest1(x);
                rst(j).fpareto{k}(:,t);
                rst(j).fpareto{k}(:,t)=f';
                rst(j).gpareto{k}(:,t)=g;
            end
            
        end
    end
    
      save(fileiter,'rst','-v7.3');    
end

%%%%%%%%%%%%%%过滤违反约束的点
for i=1:numel(filename)
    fileiter=filename{i};
    load(fileiter);
    [nrst]=length(rst);
    for  c=1:nrst
     [nc]=length(rst(c).gpareto);
     a(c).b=rst(c).gpareto{1,nc};
     a1=length(a(c).b);
     l=0;
     for j=1:a1
%          a1=length(rst(c).gpareto{1,100});
     a2=max(rst(c).gpareto{1,nc}(:,j-l));
         if max(rst(c).gpareto{1,nc}(:,j-l))>0
             rst(c).gpareto{1,nc}(:,j-l)=[];
             rst(c).ppareto{1,nc}(:,j-l)=[];
             rst(c).fpareto{1,nc}(:,j-l)=[];
             l=l+1;
         end
%          a1=length(rst(c).gpareto{1,100});
     end
             
    end
    save(fileiter,'rst','-v7.3');  
end