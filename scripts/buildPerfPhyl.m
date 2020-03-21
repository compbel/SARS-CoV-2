function AM = buildPerfPhyl(M)
n = size(M,1);
m = size(M,2);
AM = zeros(n+m+1,n+m+1);

deg = sum(M,1);
aux = [deg; 1:m]';
aux = sortrows(aux,-1);
lowest = ones(1,n);
for i=1:m
    c = aux(i,2);
    ones_c = find(M(:,c));
    if isempty(ones_c)
        continue;
    end
    v = lowest(ones_c(1));
    AM(v,c+1) = 1;
    for t=ones_c
        lowest(t) = c+1;
    end
%     
%     G = digraph(AM);
%     figure
%     plot(G,'NodeLabel',labels);
end
for t = 1:n
    AM(lowest(t),m+1+t) = 1;
end
