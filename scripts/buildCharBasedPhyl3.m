function [AM,repMuts] = buildCharBasedPhyl3(M,minorMinor)

nOrig = size(M,1);
mOrig = size(M,2);
Mobs = M;

if ~isempty(minorMinor)
    [Mperf,M,perfSeq,perfMut,nonPerfSeq,nonPerfMut] = preProcData1(M);
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
%     plot(G);
%     [comps,compsizes] = conncomp(G);
%     ind = find(comps == 1);
%     AMG = adjacency(G);
%     AMG = AMG(ind,ind);
%     figure
%     plot(graph(AMG))
    AMmut = adjacency(G);
    nonIsolMut = find(sum(AMmut,2) > 0);
    AMmut = AMmut(nonIsolMut,nonIsolMut);
    AMmutComp = ones(length(AMmut),length(AMmut)) - eye(length(AMmut),length(AMmut)) - AMmut;
    IS = maximalCliques(AMmutComp);
    VC = logical(ones(size(IS,1),size(IS,2)) - IS);
    [minsize,minset] = min(sum(VC,1));
    minsets = find(sum(VC,1) <= minsize);

    resfile = 'tree0300.mat';
    if isfile(resfile)
        load(resfile)
    else
        objs = zeros(1,length(minsets));
        repMutsOpt = [];
        AMOpt = [];
        MnperfOpt = [];
        errOpt = 1000000;
        objOpt = 1000000;
        for ms = 1:length(minsets)
%         for ms = 1:1
            minset = minsets(ms);
            repMuts = nonIsolMut(VC(:,minset))';
            mRepMut = length(repMuts);
            repSeq = unique(minorMinor(:,1));

        %     M = [0 0 0; 0 1 0; 1 0 1; 1 1 1];
        %     M = [0 0; 0 1; 1 0; 1 1];
        %     repMuts = 2;
        %     M = [1 0 0; 1 1 0; 1 0 1];
        %     repMuts = [];

%             [AM,Mnperf,err] = findTreeIPGus(M,repMuts,nonPerfMut,nonPerfSeq);
            [AM,Mnperf,err,objval] = findTreeIPGus1(M,repMuts,objOpt);
            objOpt = min(objOpt,objval);
            objs(ms) = err;
            if err < errOpt
                repMutsOpt = repMuts;
                AMOpt = AM;
                MnperfOpt = Mnperf;
                errOpt = err;
            end


        %     [AM,Mnperf] = findTreeMIP2(M,repMuts,nonPerfMut,nonPerfSeq);
        end
        AM = AMOpt;
        repMuts = repMutsOpt;
        Mnperf = MnperfOpt;
    end

    [AM,M_infer,repMuts] = extendTree(Mperf,Mnperf,Mobs,perfSeq,perfMut,nonPerfSeq,nonPerfMut,repMuts);
    
    G = digraph(AM);

    m = mOrig;
    n = nOrig;
%     drawTree(G,n,m,repMuts,[]);
    
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
        ind = find(Mtree(j,:) ~= Mobs(j,:));
        if ~isempty(ind)
            err = err + length(ind);
        end
    end
    ['end']
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