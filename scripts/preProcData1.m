function [Mperf,Mnperf,perfSeq,perfMut,nonPerfSeq,nonPerfMut] = preProcData1(M)

n = size(M,1);
m = size(M,2);

AM = buildPerfPhyl(M);
G = digraph(AM);

badSeq = zeros(1,n);
Mtree = zeros(n,m);
for j = 1:n
    path = shortestpath(G,1,m+1 + j);
    treeMut = sort(path(2:(end-1)));
    for k = treeMut
        Mtree(j,k-1) = 1;
    end
    ind = find(Mtree(j,:) ~= M(j,:));
    if ~isempty(ind)
        badSeq(j) = 1;
    end
end

perfSeq = zeros(1,n);
perfMut = zeros(1,m);
for i = 2:(m+1)
    desc = dfsearch(G,i);
    descMut = intersect(desc,1:(m+1))-1;
    descSeq = intersect(desc,(m+2):(m+1+n))-m-1;
    if sum(badSeq(descSeq),2) == 0
        nonDescSeq = setdiff(1:n,descSeq);
        if sum(sum(M(nonDescSeq,descMut),1),2) == 0
            perfSeq(descSeq) = 1;
            perfMut(descMut) = 1;
        end
    end
end

perfSeq1 = zeros(1,n);
perfMut1 = zeros(1,m);
for j = 1:n
    if perfSeq(j)
        path = shortestpath(G,1,m+1 + j);
        treeMut = sort(path(2:(end-1)))-1;
        if sum(perfMut(treeMut),2) >= length(treeMut)-0
            perfSeq1(j) = 1;
            perfMut1(treeMut) = 1;
        end
    end
end
perfMut = perfMut1;
perfSeq = perfSeq1;

noMut = find(sum(M,2) == 0);
perfSeq(noMut) = 1;
% noSeq =  find(sum(M,1) == 0);
% perfMut(noSeq) = 1;
nonPerfSeq = logical(1 - perfSeq);
nonPerfMut = logical(1 - perfMut);
Mnperf = M(nonPerfSeq,nonPerfMut);
Mperf = M(logical(perfSeq),logical(perfMut));

% idsIntern = (cellstr(num2str((0:m)')))';
% idsLeafs = (cellstr(num2str((1:n)')))';
% mutMarkers = cell(1,m+1);
% mutMarkers(:) = {'o'};
% patMarkers = cell(1,n);
% patMarkers(:) = {'s'};
% figure
% h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
% x = h.XData;
% y = h.YData;
% highlight(h,(m+2):(m+1+n),'NodeColor','r');
% highlight(h,m+1+find(perfSeq),'NodeColor','k');
% highlight(h,1+find(perfMut),'NodeColor','k');

perfSeq = find(perfSeq);
perfMut = find(perfMut);
nonPerfSeq = find(nonPerfSeq);
nonPerfMut = find(nonPerfMut);

['end'];