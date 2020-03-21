function AM = buildCharBasedPhyl2(M,minorMinor)
n = size(M,1);
m = size(M,2);

if ~isempty(minorMinor)
    M = preProcData(M);
    n = size(M,1);
    m = size(M,2);
    minorMinor = [];
    for i = 1:m
        for j = (i+1):m
            if (sum(M(:,i) + M(:,j) == 0,1) > 0) && (sum(M(:,i) - M(:,j) > 0,1) > 0) && (sum(M(:,j) - M(:,i) > 0,1) > 0) && (sum(M(:,i).*M(:,j) > 0,1) > 0)
                ind = find(M(:,i).*M(:,j) > 0);
                minorMinor = [minorMinor; [ind i*ones(length(ind),1) j*ones(length(ind),1)]];
            end
        end
     end
    
    
    edges = unique(minorMinor(:,[2 3]),'rows');
    G = graph(edges(:,1),edges(:,2));
%     figure
%     plot(G)
    AMmut = adjacency(G);
    nonIsolMut = find(sum(AMmut,2) > 0);
    AMmut = AMmut(nonIsolMut,nonIsolMut);
    AMmutComp = ones(length(AMmut),length(AMmut)) - eye(length(AMmut),length(AMmut)) - AMmut;
    IS = maximalCliques(AMmutComp);
    VC = logical(ones(size(IS,1),size(IS,2)) - IS);
    [minsize,minset] = min(sum(VC,1));

%     for s = 1:size(VC,2)
        repMuts = nonIsolMut(VC(:,minset))';
        repSeq = unique(minorMinor(:,1));
        nRepSeq = length(repSeq);
        mRepMut = length(repMuts);
        possRepMutSeq = cell(1,length(repSeq));
        for j = 1:length(repSeq)
            ind = find(minorMinor(:,1) == repSeq(j));
            mm = minorMinor(ind,[2 3]);
            possRepMutSeq{j} = intersect(unique(mm(:)),repMuts);
        end
        [possParMut{1:nRepSeq}] = ndgrid(possRepMutSeq{:});
        possParMut = reshape(cat(nRepSeq,possParMut{:}),[],nRepSeq);

        err = zeros(1,size(possParMut,1));
        badSeq = cell(1,size(possParMut,1));
        for i = 1:size(possParMut,1)
            parMutPerSeq = possParMut(i,:);
            [AM,bs,er,idsIntern] = buildExtPhyl(M,repSeq,repMuts,parMutPerSeq);
            G = digraph(AM);
            badSeq{i} = bs;
            err(i) = er;
            idsLeafs = (cellstr(num2str((1:n)')))';
            mutMarkers = cell(1,m+1+mRepMut);
            mutMarkers(:) = {'o'};
            patMarkers = cell(1,n);
            patMarkers(:) = {'s'};
            figure
            h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
            x = h.XData;
            y = h.YData;
            highlight(h,(m+mRepMut+2):(m+1+mRepMut+n),'NodeColor','r');
            highlight(h,(m+2):(m+1+mRepMut),'NodeColor','g');
            highlight(h,m+1+mRepMut+badSeq{i},'NodeColor','k');
            ['aaa'];
        end
        i = find(err == min(err));
%     end
else
        idsIntern = (cellstr(num2str((0:m)')))';
        AM = buildPerfPhyl(M);
        G = digraph(AM);
        idsLeafs = (cellstr(num2str((1:n)')))';
        mutMarkers = cell(1,m+1);
        mutMarkers(:) = {'o'};
        patMarkers = cell(1,n);
        patMarkers(:) = {'s'};
        figure
        h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
        x = h.XData;
        y = h.YData;
        highlight(h,(m+2):(m+1+n),'NodeColor','r');
        Mtree = zeros(n,m);
        for j = 1:n
            path = shortestpath(G,1,m+1 + j);
            treeMut = sort(path(2:(end-1)));
            for k = treeMut
                Mtree(j,k-1) = 1;
            end
            ind = find(Mtree(j,:) ~= M(j,:));
            if ~isempty(ind)
                b = 5;
            end
        end
end

['end']