% clear;
% fastaDir = 'sequences';
% metaDir = 'metadata';
% tree = phytreeread('nextstrain_ncov_tree.nwk');
% [A,idTree,Dist] =  getmatrix(tree);
% 
% files = dir(fullfile(fastaDir,'*.fasta'));
% allSeq = [];
% allMetadata = [];
% for f =1:size(files,1) 
%     f
%     seq_file = [fastaDir,filesep,files(f).name];
%     meta_file = [metaDir,filesep,files(f).name(1:end-5) 'pdf'];
%     seq = fastaread(seq_file);
%     name = seq.Header;
%     name = erase(name,'_');
%     name = erase(name,' ');
%     if strcmpi(name(1:8),'BetaCoV/')
%         name = name(9:end);
%     end
%     i = strfind(name,'|');
%     if i > 0
%         name = name(1:(i-1));
%     end
%     i = find(~cellfun('isempty', strfind(idTree,name)));
%     seq.Header = name;
%     allSeq = [allSeq seq];
%     if ~isempty(i)
%         seq.Header = name;
%         allSeq = [allSeq seq];
%     else
%         ['bbb'];
%     end
%     text = splitlines(extractFileText(meta_file));
%     patMeta = struct;
%     patMeta.id = name;
%     for j = 1:length(text)
%         if startsWith(text(j),"Submission Date")
%             patMeta.subdate = extractAfter(text(j),17);
%             if strlength(patMeta.date) <= 7
%                 patMeta.date = "u";
%             end
%         end
%         if startsWith(text(j),"Collection date")
%             patMeta.date = extractAfter(text(j),17);
%             if strlength(patMeta.date) <= 7
%                 patMeta.date = "u";
%         end
%         end
%         if startsWith(text(j),"Patient age")
%             patMeta.age = text(j);
%             patMeta.age = extractAfter(patMeta.age,13);
%         end
%         if startsWith(text(j),"Gender")
%             patMeta.gender = extractBetween(text(j),9,9);
%         end
%         if startsWith(text(j),"Originating lab")
%             patMeta.lab = extractAfter(text(j),16);
%         end
%         if startsWith(text(j),"Authors")
%             patMeta.authors = extractAfter(text(j),8);
%         end
%         if startsWith(text(j),"Additional host information")
%             patMeta.info = extractAfter(text(j),28);
%         end
%     end
%     allMetadata = [allMetadata patMeta];
% end
% % writetable(struct2table(allMetadata), 'gisaid.xlsx');
% 
% % subtimes = datetime(string(extractfield(allMetadata,'subdate')));
% % timesSorted = sort(subtimes);
% % timesS = unique(timesSorted);
% % nSampTime = zeros(1,length(timesS));
% % for i = 1:length(timesS)
% %     nSampTime(i) = sum(timesSorted <= timesS(i),2);
% % end
% % figure
% % plot(timesS,nSampTime);
% 
% fastawrite('nCoV.fas',allSeq);
% save('metadata.mat','allMetadata')

%% 
clear;
load('metadata.mat');
% seqF = fastaread('nCoV_align_trim_0227.fas');
% seqF = fastaread('nCoV_align_trim_0308.fas');
seqF = fastaread('nCoV_align_trim_0310.fas');
shift = 54;
refInd = 6;


seqMetadata = [];
for i = 1:length(seqF)
%     if i == 119
%         ['vvv'];
%     end
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
    else
        ['sss'];
    end
end

times = datetime(string(extractfield(seqMetadata,'date')));
seq = char(seqF.Sequence);
seqOrig = seq;
seq=fillGaps(seq);
l = size(seq,2);
n = size(seq,1);
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
% DM = l*squareform(seqpdist(seqF,'Method','p-distance'));

% remove patients: 13,74 - clustered snvs, 61 - 4 gamete, 72,76 -
% too many Ns
remove = [14,72,74,76];
seqF(remove,:) = [];
seq(remove,:) = [];
seqMetadata(remove) = [];
times(remove) = [];
ids(remove) = [];
idsShort(remove) = [];
n = n-length(remove);
ids1 = (cellstr(num2str((1:n)')))';

DM = l*squareform(pdist(seq,'hamming'));
AM = (DM <= 2) - eye(n,n);
G = graph(AM);
% plot(G,'NodeLabel',ids);

%k-step network
AMks = (DM == 0) - eye(n,n);
G = graph(AMks);
% plot(G);
[comps,compsizes] = conncomp(G);
k = 0;
while length(compsizes) > 1
    k = k + 1
    B = (DM == k);
    for i = 1:length(compsizes)
        ind = (comps == i);
        B(ind,ind) = 0;
    end
    AMks = AMks + B;
    G = graph(AMks);
%     plot(G,'NodeLabel',ids);
    [comps,compsizes] = conncomp(G);
end

% 
% figure
% plot(G,'NodeLabel',idsShort,'Layout','force');



% tree = seqlinkage(DM,'single',ids);
%% 
% % directed k-step network
% AMksdir = zeros(n,n);
% for i = 1:n
%     for j = 1:n
%         if (AMks(i,j) == 1) && (times(i) <= times(j))
%             AMksdir(i,j) = 1;
%         end
%     end
% end
% H = digraph(AMksdir);
% figure
% plot(H,'NodeLabel',idsShort,'Layout','force');

%% 
% % dates of clusters
% 
% AM1 = AMks;
% AM1(42,:) = 0;
% AM1(:,42) = 0;
% G1 = graph(AM1);
% [comps,compsizes] = conncomp(G1);
% [maxs,maxc] = max(compsizes);
% group1 = [find(comps == maxc)];
% group2 = [setdiff(1:n,group1),42];
% 
% days = min(times):max(times);
% bins1 = zeros(1,length(days));
% bins2 = zeros(1,length(days));
% for i = 1:length(days)
%     bins1(i) = sum(times(group1) == days(i),2);
%     bins2(i) = sum(times(group2) == days(i),2);
% end
% % figure
% % bar(days, [bins1; bins2])
% % G1 = graph(AMks(group1,group1));
% % figure
% % plot(G1,'NodeLabel',ids(group1));
%% 
% perfect phylogeny
% [seqAA,allMutNamesAA] = getAA(seq,shift);
% seq = seqAA;
% l = size(seq,2);

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
% fastawrite('seqvar.fas', ids,seqVar);

% getTSNE(seq,idsShort,varpos)
% getMDS(seq,idsShort)

%% 


freq = ones(1,n);
% seqVar(:,[27,32,57,72]) = [];
lVar = size(seqVar,2);
[infSiteStats,M,minorMinor] = getInfSiteStats(seqVar,freq);

% ref = seqconsensus(seqVar);
% [mint,iref] = min(times);
% ref = seqVar(iref,:);
[seqVarUn, seq_red, seq_full] = unique(seqVar, 'rows');
freqSeq = zeros(1,length(seq_red));
for i = 1:length(seq_red)
    freqSeq(i) = sum(seq_full == i,1);
end
[maxf,iref] = max(freqSeq);
ref = seqVarUn(iref,:);

M = getMutMatr(seqVar,ref);
% AMP = buildPerfPhyl(M);
% repMuts = [];
[AMP,repMuts] = buildCharBasedPhyl3(M,minorMinor);
nRepMuts = length(repMuts);
GP = digraph(AMP);

% drawTreeAA(AMP,n,lVar,repMuts,idsShort,varpos,allMutNamesAA);
% drawTree(AMP,n,lVar,repMuts,idsShort,varpos,0,[],[]);


[mutGeneLabel,mutAALabel,nonsyn] = synMut(seq,ids,shift,varpos,seq(refInd,:));
drawTree(AMP,n,lVar,repMuts,idsShort,varpos,shift,mutGeneLabel,nonsyn);

AMNet = getTransNet(AMP,GP,idsShort,times,seqMetadata,'Probhap.mat');

% mutFreq = sum(AMP(1:(lVar+nRepMuts+1),(lVar+nRepMuts+2):end),2)+1;
% unobs = find(mutFreq == 1);
% mutFreq = mutFreq/(sum(mutFreq,1));
% markSizes = 700*mutFreq;
% GPmut = digraph(AMP(1:(lVar+nRepMuts+1),1:(lVar+nRepMuts+1)));
% figure
% h = plot(GPmut,'Layout','force','NodeColor','b','MarkerSize',markSizes,'LineWidth', 2);
% highlight(h,unobs,'NodeColor','r');


% timesRel = (times - min(times))/(max(times) - min(times));
% [mint,imint] = min(times);
% [maxr,imaxt] = max(times);
% ynewL = y(imint) + timesRel.*(y(imaxt) - y(imint));
% y((lVar+2):end) = ynewL;
% figure
% h = plot(GP,'NodeLabel',[mutids ids],'Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7,'XData',x,'YData',y);
% highlight(h,(lVar+2):(lVar+1+n),'NodeColor','r');

% for i = varpos
%     nuc1 = unique(seq(group1,i));
%     nuc2 = unique(seq(group2,i));
% end
% writetable(cell2table([mutids ids]'),'perfPhylLabels.xlsx')
%% 
% [C,ia,ic] = unique(seq,'rows');
% for i = 1:length(ia)'
%     ind = find(ic == ia(i));
%     if ~isempty(ind)
%         writetable(struct2table(seqMetadata(ind)),['seq_metadata' filesep 'seq' int2str(i) '.xlsx']);
%     end
% end
%% 
% use SC method
% nCoV_param(AMP,n,lVar,times,l);

%% 
% dynamics of diversity
indNat = find(isnat(times));
times(indNat) = [];
seq(indNat,:) = [];

[seq_reg,ia,ic] = unique(seq,'rows');
iref = mode(ic);
% [mint,iref] = min(times);
ref = seq(iref,:);
timesSorted = sort(times);

% timesS = unique(timesSorted(3:end));
t1 = datetime(2020,1,1);
t2 = datetime(2020,3,10);
timesS = t1:caldays(7):t2;

divers1 = zeros(1,length(timesS));
divers2 = zeros(1,length(timesS));
divers3 = zeros(1,length(timesS));
isViol = zeros(1,length(timesS));
nVarNuc = zeros(1,length(timesS));
nNonSyn = zeros(1,length(timesS));
nSyn = zeros(1,length(timesS));
dnds = zeros(1,length(timesS));
nSeg = 10;
for i = 1:nSeg
    timesInterv(i) = timesSorted(min(length(timesSorted),round(i*length(timesSorted)/nSeg)));
    ind = find(times <= timesInterv(i));
% for i = 1:length(timesS)
%     ind = find(times <= timesS(i));
    distRef = l*pdist2(seq(ind,:),ref,'hamming');
    divers1(i) = mean(distRef);
    seq1 = unique(seq(ind,:),'rows');
    distRef = l*pdist2(seq1,ref,'hamming');
    divers2(i) = mean(distRef);
    distAv = l*pdist(seq1,'hamming');
    divers3(i) = mean(distAv);
    [infSiteStats,M,minorMinor] = getInfSiteStats(seqVar(ind,:),freq(ind));
    isViol(i) = size(minorMinor,1);
    
    indVar = zeros(1,l);
    for j = 1:l
        un = unique(seq(ind,j));
        if length(un) > 1
            nVarNuc(i) = nVarNuc(i)+1;
            indVar(j) = 1;
        end
    end
    varpos = find(indVar == 1);
    [mutGeneLabel,nonsyn] = synMut(seq1,ids,shift,varpos);
    nNonSyn(i) = length(nonsyn);
    nSyn(i) = length(varpos) - nNonSyn(i);
    dnds(i) = nNonSyn(i)/nSyn(i);
end
% timesS = timesInterv;
figure
plot(timesS,divers2)
figure
plot(timesS,divers3)
figure
plot(timesS,isViol)
figure
plot(timesS,nVarNuc)
figure
plot(timesS,nNonSyn)
figure
plot(timesS,nSyn)
figure
plot(timesS,dnds)
%% 
% % centrality
% deg = sum(AMks,2);
% centr = centrality(G,'betweenness');
% ind = find(centr > 50);
% md = seqMetadata(ind);
% 
% tree = phytreeread('nextstrain_ncov_tree.nwk');
% % bintree = phytree2bintree(tree,n,l,ids);
% % AM_tree = treeAMdirect(bintree);
% % G = digraph(AM_tree);
% % figure
% % plot(G);