function [x,t,var] = tree2SolIP(AM,n,m)
x = zeros(n,m);
t = zeros(m,m);

varInd = [n*m m*(m-1)/2];
nVar = sum(varInd,2);
var = zeros(1,nVar);

G = digraph(AM);

for i = 1:n
    path = shortestpath(G,1,1+m+i);
    for j=path(2:end-1)
        x(i,j-1)=1;
    end
end

for j1=1:m
    for j2=(j1+1):m
        if sum(x(:,j1).*x(:,j2),1) > 0
            if sum(x(:,j1) - x(:,j2),1) > 0
                t(j1,j2) = 1;
            else
                t(j2,j1) = 1;
            end
        end
    end
end

count = 1;
for i = 1:n
    for j = 1:m
        var(count) = x(i,j);
        count = count + 1;
    end
end

for j1 = 1:m
    for j2=(j1+1):m
        var(count) = t(j1,j2);
        count = count+1;
    end
end