clear;
% loads sequence alignment and metadata
load('metadata.mat');
seqF = fastaread('nCoV_align_trim_0310.fas');

% parse sequence IDs from metadata
seqMetadata = [];
for i = 1:length(seqF)
    seqF(i).Header = erase(seqF(i).Header,'_');
    seqF(i).Header = erase(seqF(i).Header,' ');
    found = false;
    for j = 1:length(allMetadata)
        md = allMetadata(j);
        mdName = md.id;
        if strcmp(mdName(1:9),' BetaCoV/')
            mdName = mdName(10:end);
        end
        if strcmp(seqF(i).Header,mdName)
            found = true;
            break;
        end
    end
    if found
        seqMetadata = [seqMetadata md];
    end
end

% extract sampling times
times = datetime(string(extractfield(seqMetadata,'date')));
seq = char(seqF.Sequence);
seqOrig = seq;
seq=fillGaps(seq);
l = size(seq,2);
n = size(seq,1);
%parse names of geographical locations from sequence IDs
ids = extractfield(seqF,'Header');
idsShort = cell(1,n);
for i = 1:length(ids)
    ids{i} = [ids{i}(1:end-4) char(seqMetadata(i).date) '/'];
    if strcmpi(ids{i}(1:4),'hCoV')
        ids{i} = ids{i}(9:end);
    end
    idsTok = split(ids{i},["/","-"]);
    idsShort{i} = idsTok{1};
    if strcmp(idsShort{i},'USA')
        nextTok = idsTok{2};
        if startsWith(nextTok,'Cruise')
             idsShort{i} = [idsShort{i} '-' 'Cruise'];
        else      
            idsShort{i} = [idsShort{i} '-' nextTok(1:2)];
        end
    end
    ids{i} = strrep(ids{i},'_','-');
    idsShort{i} = strrep(idsShort{i},'_','-');
end


% number of bootstraps
nBS = 1000;

%extract variable nucleotide positions
ind = ones(1,l);
varNuc = [];
for i = 1:l
    un = unique(seq(:,i));
    if length(un) == 1
        ind(i) = 0;
    else
        varNuc = [varNuc length(un)];
    end
end
varpos = find(ind == 1);
seqVar = seq(:,varpos);
seqOrigVar = seqOrig(:,varpos);
freq = ones(1,n);
lVar = size(seqVar,2);

% bootstrapping
for bs = 1:nBS  
    ind = sort(datasample(1:lVar,lVar));
    seqVarB = seqVar(:,ind);

    % identifies 4-gamete rule violations
    [infSiteStats,M,minorMinor] = getInfSiteStats(seqVarB,freq);

    % identifies the reference as the most frequent sequence
    [seqVarUn, seq_red, seq_full] = unique(seqVarB, 'rows');
    freqSeq = zeros(1,length(seq_red));
    for i = 1:length(seq_red)
        freqSeq(i) = sum(seq_full == i,1);
    end
    [maxf,iref] = max(freqSeq);
    ref = seqVarUn(iref,:);

    % constructs mutation matrix
    M = getMutMatr(seqVarB,ref);
    % construct character-based phylogeny
    [AMP,repMuts] = buildCharBasedPhyl3(M,minorMinor);
    GP = digraph(AMP);
    AMNet = getTransNet0(AMP,GP,idsShort,times,seqMetadata);
    % saves the results of the given bootstrap to a file
    save(['bootstraping_mut' filesep 'bootstrap' int2str(bs) '.mat']);
end
%% 
clear;
load(['bootstraping_mut' filesep 'bootstrap1.mat']);
n = length(idsShort);
nBoot = nBS;
AMBoot = zeros(n,n);
for bs = 1:nBoot  
    bs
    % loads results of a given bootstrap 
    load(['bootstraping_mut' filesep 'bootstrap' int2str(bs) '.mat']);
    G = digraph(AMP);
    m = size(AMP,1) - n - 1;
    parLeafs = zeros(1,n);
    for i = 1:n
        parLeafs(i) = find(AMP(:,m+1+i));
    end
    haps = unique(parLeafs);
    nHaps = length(haps);

    % constructs potential transmission links between sequences corresponding to the given bootstrap
    H = transclosure(G);
    AMcl = adjacency(H);
    AMcl = AMcl(haps,haps);
    H = transreduction(digraph(AMcl));
    AMNet = adjacency(H);
    % counts how many times an edge was generated during the bootstrapping
    for i = 1:n
        for j = 1:n
            if i == j
                continue;
            end
            u = find(haps == parLeafs(i));
            v = find(haps == parLeafs(j));
            if (AMNet(u,v) == 1) || (u == v)
                AMBoot(i,j) = AMBoot(i,j) + 1;
            end
        end
    end
end

% normalize edge counts
AMBoot = AMBoot/nBS;

%% 
% loads haplotypes (the sets of identical sequences) and their attributes
load('haplotypes0310.mat');
nHaps = length(hapStruct);
% calculates probabilities of edges between haplotypes
Probhap = zeros(nHaps,nHaps);
for i = 1:nHaps
    for j = 1:nHaps
        if i == j
            continue;
        end
        Probhap(i,j) = max(max(AMBoot(hapStruct(i).Sequences,hapStruct(j).Sequences)));
    end
end

% builds transmission network
AMNet = (Probhap > 0.75);
G = digraph(AMNet);
[comps,compsizes] = conncomp(G,'Type','weak');
[maxsize,maxcomp] = max(compsizes);
maxcomp = find(comps == maxcomp);
while maxsize < nHaps
    roots = [];
    for i = 1:length(compsizes)
        if compsizes(i) == maxsize
            continue;
        end
        comp = find(comps == i);
        indeg = sum(AMNet(comp,comp),1);
        r = comp(indeg == 0);
        roots = [roots r];
    end
    P = Probhap(maxcomp,roots);
    [u,v] = find(P == max(max(P)));
    for j = 1:length(u)
        AMNet(maxcomp(u(j)),roots(v(j))) = 1;
    end
    G = digraph(AMNet);
    [comps,compsizes] = conncomp(G,'Type','weak');
    [maxsize,maxcomp] = max(compsizes);
    maxcomp = find(comps == maxcomp);
end

% plots transmission network
H = digraph(AMNet.*Probhap);
figure
h = plot(H,'NodeLabel',hapsID,'Layout','force','NodeColor','b','MarkerSize',7,'EdgeLabel',H.Edges.Weight);
h.NodeFontSize = 7;

% saves node and edge attributes to be used in Gephi
hapStruct = rmfield(hapStruct,'Sequences');
hapStruct = rmfield(hapStruct,'Metadata');
writetable(struct2table(hapStruct), 'transNetAttrNBS.csv');
[s,t] = find(AMNet);
for i = 1:length(s)
    edge(i).Source = s(i);
    edge(i).Target = t(i);
    edge(i).Weight = Probhap(s(i),t(i));
end
writetable(struct2table(edge), 'transNetAttrEBS.csv');
