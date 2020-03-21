function AM = buildCharBasedPhyl(M,minorMinor)
n = size(M,1);
m = size(M,2);
AM = zeros(n+m+1,n+m+1);

% labels = cell(1,n+m+1);
% for i = 0:m
%     labels{i+1} = ['m_' int2str(i)];
% end
% for i = 1:n
%     labels{m+1+i} = ['c_' int2str(i)];
% end

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

corrProf = zeros(1,n);
falsePos = cell(1,n);
falseNeg = cell(1,n);
G = digraph(AM);
for i = 1:n
    path = shortestpath(G,1,m+1+i);
    treeMut = sort(path(2:(end-1))-1);
    trueMut = sort(find(M(i,:) > 0));
    if isequal(treeMut,trueMut)
        corrProf(i) = 1;
    end
    falsePos{i} = setdiff(treeMut,trueMut);
    falseNeg{i} = setdiff(trueMut,treeMut);
end
falsePosSeq = find(~cellfun(@isempty,falsePos));
falseNegSeq = find(~cellfun(@isempty,falseNeg));

idsIntern = (cellstr(num2str((0:m)')))';
idsLeafs = (cellstr(num2str((1:n)')))';
AMext = AM;
mold = m;
bflab = distances(G,1);
nodeMutMap = 1:m;
for i = falseNegSeq
    row = find(minorMinor(:,1) == i);
    if length(falseNeg{i}) == 1
        falseNegMut = falseNeg{i};
        possRepMuts = [];
        for j = row'
            if minorMinor(j,2) == falseNegMut
                possRepMuts = [possRepMuts minorMinor(j,3)];
            else
                possRepMuts = [possRepMuts minorMinor(j,2)];
            end
        end
        [maxL,imaxL] = max(bflab(possRepMuts+1));
        repMut = possRepMuts(imaxL);
        nodeMutMap = [nodeMutMap repMut];
        path = shortestpath(G,repMut+1,m+1+i);
        v = path(2);
        if v > m+1
            v = v+1;
        end
        AMext1 = zeros(n+m+2,n+m+2);
        AMext1(1:(m+1),1:(m+1)) = AMext(1:(m+1),1:(m+1));
        AMext1(1:(m+1),(m+3):end) = AMext(1:(m+1),(m+2):end);
        AMext1(repMut+1,v) = 0;
        AMext1(falseNegMut+1,m+2) = 1;
        AMext1(m+2,v) = 1;
        m = m + 1;
        AMext = AMext1;
        idsIntern{end+1} = [int2str(repMut) 'r'];
        G = digraph(AMext);
        bflab = distances(G,1);
        
    end
end

mutMarkers = cell(1,m+1);
mutMarkers(:) = {'o'};
patMarkers = cell(1,n);
patMarkers(:) = {'s'};
figure
h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
x = h.XData;
y = h.YData;
highlight(h,(m+2):(m+1+n),'NodeColor','r');
highlight(h,(mold+2):(m+1),'NodeColor','g');

corrProf = zeros(1,n);
falsePos = cell(1,n);
falseNeg = cell(1,n);
G = digraph(AMext);
for i = 1:n
    path = shortestpath(G,1,m+1+i);
    treeMut = sort(unique(nodeMutMap(path(2:(end-1))-1)));
    trueMut = sort(find(M(i,:) > 0));
    if isequal(treeMut,trueMut)
        corrProf(i) = 1;
    end
    falsePos{i} = setdiff(treeMut,trueMut);
    falseNeg{i} = setdiff(trueMut,treeMut);
end
falsePosSeq = find(~cellfun(@isempty,falsePos));
falseNegSeq = find(~cellfun(@isempty,falseNeg));

for i = falsePosSeq
    falsePosMut = falsePos{i};
end

['end']