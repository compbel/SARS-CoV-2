function AM = sol2TreeIP(x)
n = size(x,1);
m = size(x,2);
AM = zeros(n+m+1,n+m+1);

% labels = cell(1,n+m+1);
% for i = 0:m
%     labels{i+1} = ['m_' int2str(i)];
% end
% for i = 1:n
%     labels{m+1+i} = ['c_' int2str(i)];
% end

deg = sum(x,1);
aux = [deg; 1:m]';
aux = sortrows(aux,-1);
lowest = ones(1,n);
for i=1:m
    c = aux(i,2);
    ones_c = find(x(:,c));
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