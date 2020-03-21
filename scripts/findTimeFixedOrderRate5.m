function [times,likel] = findTimeFixedOrderRate5(stree,AM,mutOrders,rates,minRate,maxRate)
% rate-related constraints on differences between birth and last child
% birth; no sampling time
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
leafs = find(deg == 0);
cIntern = (deg > 0); 
intern = find(deg > 0);
nLeafs = length(leafs);
nIntern = length(intern);


A = zeros(m-1+2*nIntern + nLeafs,m+1);
rhs = zeros(m-1+2*nIntern + nLeafs,1);

iEq = 0;
for i = 1:m
    if ~isempty(stree(i).children)
        for j = 2:length(mutOrders{i})
            iEq = iEq + 1;
            A(iEq,mutOrders{i}(j-1)) = 1;
            A(iEq,mutOrders{i}(j)) = -1;
        end
    end
end

for i = intern
    iEq = iEq + 1;
    j = mutOrders{i}(end);
    A(iEq,i) = 1;
    A(iEq,j) = -1;
    rhs(iEq) = -deg(i)/maxRate;
    
    iEq = iEq + 1;
    A(iEq,i) = -1;
    A(iEq,j) = 1;
    rhs(iEq) = deg(i)/minRate;
end

for i = leafs
    iEq = iEq + 1;
    A(iEq,i) = 1;
    A(iEq,m+1) = -1;
end


nEq = 0;
for i = intern
    if ismember(mutOrders{i}(end),leafs)
        nEq = nEq + 1;
    end
end
Aeq = zeros(nEq,m+1);
rhseq = zeros(nEq,1);
iEq = 1;
for i = intern
    if ismember(mutOrders{i}(end),leafs)
        Aeq(iEq, mutOrders{i}(end)) = 1;
        Aeq(iEq, m+1) = -1;
        iEq = iEq + 1;
    end
end


lb = zeros(1,m+1);
ub = (m-1)/minRate*ones(1,m+1);

% lb = deg/maxRate;
% ub = deg/minRate;
% ub(leafs) = (m-1)/minRate;
ub(1) = 0;

theta0 = 0.5*(maxRate + minRate);
currTime = zeros(1,m);
AMOrders = zeros(m,m);
    for i = 1:m
        for j = 2:length(mutOrders{i})
            AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
        end
    end
G = digraph(AMOrders);
% G = digraph(AM);
% v = bfsearch(G,1);
% for i = 1:m
%     currTime(v(i)) = (i-1)/theta0;
% end

tStep = 10;
or = toposort(G);
step = 0;
for i=or
    currTime(i) = step*tStep;
    step = step +1;
end


currTime(m+1) = max(currTime(leafs));

% currTheta = extractfield(stree,'rate');
% currTime = extractfield(stree,'time');
% currTime(1) = 0.00001;

% options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',false);
% options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',15000);
% l = @(t)likelihoodRatesOrders1(t,rates,mutOrders);
% [times,likel] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);


[currTime,EMOrders,A1,rhs1] = getInitTimeConsMatr(mutOrders,stree,minRate,maxRate,10000,m,leafs);

    mu = [rates.*cIntern, 0];
    lambda = zeros(1,m+1);
    for i = intern
        lambda(mutOrders{i}(end)) = rates(i);
    end
%     cvx_expert true
    cvx_begin quiet
%     cvx_begin
        variable t(m+1)
        minimize -mu*t + lambda*t - sum(log(EMOrders*t))
        A*t <= rhs;
        Aeq*t == rhseq;
%         Abd*t <= rhsbd
        lb' <= t <= ub'
    cvx_end
    likel= likelihoodRatesOrders1(t,rates,mutOrders);
    times = t';