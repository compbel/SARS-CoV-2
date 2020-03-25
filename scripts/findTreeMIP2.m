function [AM,M_infer_full] = findTreeMIP2(M,repMut,nonPerfMut,nonPerfSeq)
eps = 0.01;

M_full = M;
[M, mut_red, mut_full] = unique(M_full', 'rows');
M = M';
repMut_full = repMut;
repMut = mut_full(repMut)';
[M, seq_red, seq_full] = unique(M, 'rows');

n = size(M,1);
mOrig = size(M,2);
mRepMut = length(repMut);
m = mOrig+mRepMut;
allMuts = [1:mOrig repMut];


% [row, col] = find(M == 0);
% E1 = [row, col];
% [row, col] = find(M == 1);
% E2 = [row, col];

% variables: 1)x_{i,j} 2)u_i 3)v_i 2)t_{i,j}  
varInd = [(n+m)*(n+m) n+m n+m (n+m)*(n+m-1)/2];
nVar = sum(varInd,2);
nEq = 2*(n+m)*(n+m-1) + (n+m)*(n+m-1) + (n+m)*(n+m-1) + (n+m) + m*(m-1)/2

indx = 0;
indu = varInd(1);
indv = varInd(1) + varInd(2);
indt = varInd(1) + varInd(2) + varInd(3);

%objective
f = zeros(nVar,1);
% for j = (n+1):(n+m)
%     for i = 1:n
%             j1 = allMuts(j - n);
%             if M(i,j1) == 0
%                 f((j-1)*(n+m)+i) = 1;
%             else
%                 f((j-1)*(n+m)+i) = -1;
%             end
%     end
% end

for i = 1:n
    for j = 1:mOrig
        ind = find(allMuts == j);
        if M(i,j) == 0
            for j1 = ind
                j2 = n+j1;
                f((j2-1)*(n+m)+i) = 1;
            end
        else
            for j1 = ind
                j2 = n+j1;
                f((j2-1)*(n+m)+i) = -1;
            end
        end
    end
end

A = zeros(nEq,nVar);
rhs = zeros(nEq,1);
eq = 1;

% u_i - u_j \leq 1-x_{i,j}; v_j - v_i \leq 1 - x_{i,j}
for i = 1:(n+m)
    for j = 1:(n+m)
        if i == j
            continue;
        end
        A(eq,indu + i) = 1;
        A(eq,indu + j) = -1;
        A(eq,indx + (n+m)*(i-1) + j) = 1;
        rhs(eq) = 1;
        eq = eq + 1;
        A(eq,indv + i) = -1;
        A(eq,indv + j) = 1;
        A(eq,indx + (n+m)*(i-1) + j) = 1;
        rhs(eq) = 1;
        eq = eq + 1;
    end
end

% -t_{i,j} \leq u_i - u_j \leq 1 - t_{i,j)
it = 1;
for i = 1:(n+m)
    for j = (i+1):(n+m)
        A(eq,indt+it) = -1;
        A(eq,indu + i) = -1;
        A(eq,indu + j) = 1;
        rhs(eq) = 0;
        eq = eq + 1;
        A(eq,indt+it) = 1;
        A(eq,indu + i) = 1;
        A(eq,indu + j) = -1;
        rhs(eq) = 1;
        eq = eq + 1;
        it = it+1;
    end
end

% v_i - u_j \leq x_{i,j} + x_{j,i} + 1 - t_{i,j}; v_j - u_i \leq x_{i,j} +
% x_{j,i} + t_{i,j}

it = 1;
for i = 1:(n+m)
    for j = (i+1):(n+m)
        A(eq,indv + i) = 1;
        A(eq,indu + j) = -1;
        A(eq,indx + (i-1)*(n+m) + j) = -1;
        A(eq,indx + (j-1)*(n+m) + i) = -1;
        A(eq,indt + it) = 1;
        rhs(eq) = 1;
        eq = eq + 1;
        A(eq,indv + j) = 1;
        A(eq,indu + i) = -1;
        A(eq,indx + (i-1)*(n+m) + j) = -1;
        A(eq,indx + (j-1)*(n+m) + i) = -1;
        A(eq,indt + it) = -1;
        rhs(eq) = 0;
        eq = eq + 1;
        it = it+1;
    end
end

% u_i \leq v_i
for i = 1:(n+m)
    A(eq,indu + i) = 1;
    A(eq,indv + i) = -1;
    rhs(eq) = -eps;
    eq = eq + 1;
end


% x_{i,j} + x_{j,i} \leq 1;
for i = (n+1):(n+m)
    for j = (i+1):(n+m)
        A(eq,indx + (i-1)*(n+m) + j) = 1;
        A(eq,indx + (j-1)*(n+m) + i) = 1;
        rhs(eq) = 1;
        eq = eq + 1;
    end
end

%x_{i,j} = 0 for i = 1:n

lb = zeros(nVar,1);
ub = ones(nVar,1);
for i = 1:n
    for j = 1:(n+m)
         ub(indx + (i-1)*(n+m) + j) = 0;
    end
end

for j = 1:m
    ind = find(allMuts == j);
    if length(ind) == 2
        u = n+ind(1);
        v = n+ind(2);
        ub(indx + (u-1)*(n+m) + v) = 0;
        ub(indx + (v-1)*(n+m) + u) = 0;
    end
end

varType = char(ones(1,nVar)*'B');
for i = 1:(n+m)
    varType(indu + i) = 'C';
    varType(indv + i) = 'C';
end


clear model;
model.A = sparse(A);
model.obj = f';
model.rhs = rhs;
model.vtype = varType;
model.sense = '<';
model.modelsense = 'min';
model.lb = lb;
model.ub = ub;
% model.bestobjstop = ld + 0.1;
params.outputflag = 1;
params.timelimit = 1000;
% params.MIPFocus = 2;
result = gurobi(model, params);

[M_infer,x,u,v,t] = getMIPsolLPInterv(result.x,n,m);

% load('test1_red.mat');


%     m = size(M,2);
%     AM = buildPerfPhyl(M_infer);
%     G = digraph(AM);
% %     idsIntern = (cellstr(num2str((0:m)')))';
% %     for j = 1:mRepMut
% %         idsIntern{end+1} = [int2str(repMut(j)) 'r'];
% %     end
% %     idsLeafs = (cellstr(num2str((1:n)')))';
%     idsIntern = (cellstr(num2str([0 mut_red']')))';
%     for j = 1:mRepMut
%         idsIntern{end+1} = [int2str(repMut_full(j)) 'r'];
%     end
%     idsLeafs = (cellstr(num2str((1:n)')))';
%     mutMarkers = cell(1,m+1+mRepMut);
%     mutMarkers(:) = {'o'};
%     patMarkers = cell(1,n);
%     patMarkers(:) = {'s'};
%     figure
%     h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
%     highlight(h,(m+mRepMut+2):(m+1+mRepMut+n),'NodeColor','r');
%     highlight(h,(m+2):(m+1+mRepMut),'NodeColor','g');
    
M_infer_full = zeros(size(M_full,1),size(M_full,2)+mRepMut);
M_infer = M_infer(seq_full,:);
for j = 1:size(M_full,2)
    M_infer_full(:,j) = M_infer(:,mut_full(j));
end
for j = 1:mRepMut
    M_infer_full(:,size(M_full,2) + j) = M_infer(:,size(M,2)+j);
end

m = size(M_full,2);
n = size(M_full,1);
AM = buildPerfPhyl(M_infer_full);
G = digraph(AM);
Mtree = zeros(n,m);
err = 0;
for j = 1:n
    path = shortestpath(G,1,m+1+mRepMut + j);
    treeMut = sort(path(2:(end-1)));
    for k = treeMut
        if k <= m+1
            Mtree(j,k-1) = 1;
        else
            Mtree(j,repMut_full(k-m-1)) = 1;
        end
    end
    ind = find(Mtree(j,:) ~= M_full(j,:));
    if ~isempty(ind)
        err = err + length(ind);
    end
end
    

    
['end']
end
    

