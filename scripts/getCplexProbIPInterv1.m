function [f,A,rhs,lb,ub,varType] = getCplexProbIPInterv1(M,repMuts)
eps = 0.01;

n = size(M,1);
mOrig = size(M,2);
mRepMuts = length(repMuts);
m = mOrig+mRepMuts;
allMuts = [1:mOrig repMuts];

% [row, col] = find(M == 0);
% E1 = [row, col];
% [row, col] = find(M == 1);
% E2 = [row, col];

% variables: 1)x_{i,j} 2)u_i 3)v_i 2)t_{i,j}  
varInd = [(n+m)*(n+m) n+m n+m (n+m)*(n+m-1)/2];
nVar = sum(varInd,2);
nEq = 2*(n+m)*(n+m-1) + (n+m)*(n+m-1) + (n+m)*(n+m-1) + (n+m) + m*(m-1)/2 + n*(n+m)

indx = 0;
indu = varInd(1);
indv = varInd(1) + varInd(2);
indt = varInd(1) + varInd(2) + varInd(3);

%objective
f = zeros(nVar,1);
for j = (n+1):(n+m)
    for i = 1:n
            j1 = allMuts(j - n);
            if M(i,j1) == 0
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
for i = 1:n
    for j = 1:(n+m)
         A(eq,indx + (i-1)*(n+m) + j) = 1;
         rhs(eq) = 0;
         eq = eq + 1;
    end
end

lb = zeros(nVar,1);
ub = ones(nVar,1);
varType = [];

varType = char(ones(1,nVar)*'B');
for i = 1:(n+m)
    varType(indu + i) = 'C';
    varType(indv + i) = 'C';
end


['end']
end
    

