function M_infer = findTreeMIP(M,repMuts)
eps = 0.01;

n = size(M,1);
mOrig = size(M,2);
Mext = [ones(n,1) M];
repMutsExt = repMuts + 1;
mRepMut = length(repMuts);
m = mOrig+mRepMut+1;
allMutsExt = [1:(mOrig+1) repMutsExt];


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
for j = (n+2):(n+m)
    for i = 1:n
            j1 = allMutsExt(j - n);
            if Mext(i,j1) == 0
                f((j-1)*(n+m)+i) = 1;
            else
                f((j-1)*(n+m)+i) = -1;
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
params.timelimit = 1200;
params.MIPFocus = 2;
result = gurobi(model, params);

[M_infer_ext,x,u,v,t] = getMIPsolLPInterv(result.x,n,m);
m = mOrig;
M_infer = M_infer_ext(:,2:end);
AM = buildPerfPhyl(M_infer);
G = digraph(AM);
Mtree = zeros(n,m);
err = 0;
for j = 1:n
    path = shortestpath(G,1,m+1+mRepMut + j);
    treeMut = sort(path(2:(end-1)));
    for k = treeMut
        if k <= mOrig+1
            Mtree(j,k-1) = 1;
        else
            Mtree(j,repMuts(k-m-1)) = 1;
        end
    end
    ind = find(Mtree(j,:) ~= M(j,:));
    if ~isempty(ind)
        err = err + length(ind);
    end
end
obj = 0;
    for u = 1:n
        for v = 1:m
            if M(u,v) == 0
                obj = obj + Mtree(u,v);
            else
                obj = obj - Mtree(u,v);
            end
        end
    end

    
    idsIntern = (cellstr(num2str((0:m)')))';
    for j = 1:mRepMut
        idsIntern{end+1} = [int2str(repMuts(j)) 'r'];
    end
    idsLeafs = (cellstr(num2str((1:n)')))';
    mutMarkers = cell(1,m+1+mRepMut);
    mutMarkers(:) = {'o'};
    patMarkers = cell(1,n);
    patMarkers(:) = {'s'};
    figure
    h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
    highlight(h,(m+mRepMut+2):(m+1+mRepMut+n),'NodeColor','r');
    highlight(h,(m+2):(m+1+mRepMut),'NodeColor','g');
    
['end']
end
    

