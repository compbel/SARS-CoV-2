% clear;
% load('metadata.mat');
% seqF = fastaread('nCoV_align_trim_0310.fas');
% shift = 54;
% refInd = 6;
% 
% 
% seqMetadata = [];
% for i = 1:length(seqF)
% %     if i == 119
% %         ['vvv'];
% %     end
%     seqF(i).Header = erase(seqF(i).Header,'_');
%     seqF(i).Header = erase(seqF(i).Header,' ');
%     found = false;
%     for j = 1:length(allMetadata)
%         md = allMetadata(j);
%         mdName = md.id;
%         if strcmp(mdName(1:9),' BetaCoV/')
%             mdName = mdName(10:end);
%         end
%         if strcmp(seqF(i).Header,mdName)
%             found = true;
%             break;
%         end
%     end
%     if found
%         seqMetadata = [seqMetadata md];
%     else
%         ['sss'];
%     end
% end
% 
% times = datetime(string(extractfield(seqMetadata,'date')));
% seq = char(seqF.Sequence);
% seqOrig = seq;
% seq=fillGaps(seq);
% l = size(seq,2);
% n = size(seq,1);
% ids = extractfield(seqF,'Header');
% idsShort = cell(1,n);
% for i = 1:length(ids)
%     ids{i} = [ids{i}(1:end-4) char(seqMetadata(i).date) '/'];
%     if strcmpi(ids{i}(1:4),'hCoV')
%         ids{i} = ids{i}(9:end);
%     end
%     idsTok = split(ids{i},["/","-"]);
%     idsShort{i} = idsTok{1};
%     if strcmp(idsShort{i},'USA')
%         nextTok = idsTok{2};
%         if startsWith(nextTok,'Cruise')
%              idsShort{i} = [idsShort{i} '-' 'Cruise'];
%         else      
%             idsShort{i} = [idsShort{i} '-' nextTok(1:2)];
%         end
%     end
%     ids{i} = strrep(ids{i},'_','-');
%     idsShort{i} = strrep(idsShort{i},'_','-');
% end
% 
% remove = [14,72,74,76];
% seqF(remove,:) = [];
% seq(remove,:) = [];
% seqMetadata(remove) = [];
% times(remove) = [];
% ids(remove) = [];
% idsShort(remove) = [];
% n = n-length(remove);
% ids1 = (cellstr(num2str((1:n)')))';
% 
% nBS = 500;
% 
% ind = ones(1,l);
% varNuc = [];
% for i = 1:l
%     un = unique(seq(:,i));
%     if length(un) == 1
%         ind(i) = 0;
%     else
%         varNuc = [varNuc length(un)];
%     end
% end
% varpos = find(ind == 1);
% seqVar = seq(:,varpos);
% seqOrigVar = seqOrig(:,varpos);
% freq = ones(1,n);
% lVar = size(seqVar,2);
%     
% for bs = 258:300  
%     ind = sort(datasample(1:lVar,lVar));
%     seqVarB = seqVar(:,ind);
% %     ind = sort(datasample(1:n,n));
% %     seqVarB = seqVar(ind,:);
%     
%     [infSiteStats,M,minorMinor] = getInfSiteStats(seqVarB,freq);
% 
%     [seqVarUn, seq_red, seq_full] = unique(seqVarB, 'rows');
%     freqSeq = zeros(1,length(seq_red));
%     for i = 1:length(seq_red)
%         freqSeq(i) = sum(seq_full == i,1);
%     end
%     [maxf,iref] = max(freqSeq);
%     ref = seqVarUn(iref,:);
% 
%     M = getMutMatr(seqVarB,ref);
%     % AMP = buildPerfPhyl(M);
%     % repMuts = [];
%     [AMP,repMuts] = buildCharBasedPhyl3(M,minorMinor);
%     GP = digraph(AMP);
%     AMNet = getTransNet0(AMP,GP,idsShort,times,seqMetadata);
%     save(['bootstraping_mut' filesep 'bootstrap' int2str(bs) '.mat']);
% end
%% 
clear;
load(['bootstraping_mut' filesep 'bootstrap1.mat']);
n = length(idsShort);
AMpatBoot = zeros(n,n);
nBS1 = 300;
for bs = 1:nBS1  
    bs
    load(['bootstraping_mut' filesep 'bootstrap' int2str(bs) '.mat']);
    G = digraph(AMP);
    m = size(AMP,1) - n - 1;
    parLeafs = zeros(1,n);
    for i = 1:n
        parLeafs(i) = find(AMP(:,m+1+i));
    end
    haps = unique(parLeafs);
    nHaps = length(haps);

    H = transclosure(G);
    AMcl = adjacency(H);
    AMcl = AMcl(haps,haps);
    H = transreduction(digraph(AMcl));
    AMNet = adjacency(H);
    for i = 1:n
        for j = 1:n
            if i == j
                continue;
            end
            u = find(haps == parLeafs(i));
            v = find(haps == parLeafs(j));
            if (AMNet(u,v) == 1) || (u == v)
                AMpatBoot(i,j) = AMpatBoot(i,j) + 1;
            end
        end
    end
end

AMpatBoot = AMpatBoot/nBS1;

%% 


%imported from getTransNet
load('haplotypes0310.mat');
nHaps = length(hapStruct);
Probhap = zeros(nHaps,nHaps);
for i = 1:nHaps
    for j = 1:nHaps
        if i == j
            continue;
        end
        Probhap(i,j) = max(max(AMpatBoot(hapStruct(i).Sequences,hapStruct(j).Sequences)));
    end
end

AMNet = (Probhap > 0.75);
G = digraph(AMNet);
% figure
% plot(G,'Layout','force','NodeLabel',hapsID,'NodeColor','b','MarkerSize',7)
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
%     figure
%     plot(G,'Layout','force','NodeLabel',hapsID,'NodeColor','b','MarkerSize',7)
    [comps,compsizes] = conncomp(G,'Type','weak');
    [maxsize,maxcomp] = max(compsizes);
    maxcomp = find(comps == maxcomp);
end
H = digraph(AMNet.*Probhap);
figure
h = plot(H,'NodeLabel',hapsID,'Layout','force','NodeColor','b','MarkerSize',7,'EdgeLabel',H.Edges.Weight);
% h = plot(H,'NodeLabel',hapsID,'Layout','force','NodeColor','b','MarkerSize',7);
h.NodeFontSize = 7;
ind = find(contains(hapsID,'Singapore'));
highlight(h,ind,'NodeColor','r');

intraProb = [];
for i = ind
    for j = ind
        if AMNet(i,j) == 1
            intraProb = [intraProb Probhap(i,j)];
        end
    end
end
        

% ind1 = find(contains(hapsID,'USA-WA'));
% ind = setdiff(ind,ind1);

% dist = distances(H,1);
% mean(dist(ind))
% std(dist(ind))

load('countries.mat');

AMNetUndir = (AMNet + AMNet') > 0;
deg = sum(AMNetUndir,2);
dd = histc(deg, unique(deg));
figure
loglog(dd,'LineWidth',1);
xlabel('degree')
ylabel('count')

pd = fitdist(dd,'GeneralizedPareto')
y = pdf(pd,1:length(dd));
figure
plot(1:length(dd),dd/sum(dd,1),1:length(dd),y,'LineWidth',1);
% plot(1:length(dd),y,'LineWidth',1);


k = find(dd == 1,1);
coeffs = polyfit(1:k,log(dd(1:k))', 1)
figure
plot(1:k,log(dd(1:k))');


hapsTimes = [];
for i = 1:nHaps
    t = min(times(parLeafs == haps(i)));
    hapsTimes = [hapsTimes t];
end
Dist = distances(H);
D1 = [];
D2 = [];
hapsTimesNum = days(hapsTimes - datetime(2019,12,1));
for i = 1:nHaps
    v = find(~isinf(Dist(i,:)));
    D1 = [D1 Dist(i,v)];
    D2 = [D2  (hapsTimesNum(v)-(hapsTimesNum(i)+2)*ones(1,length(v)))];
end
ind = find(~isnan(D2));
D1 = D1(ind);
D2 = D2(ind);
[rho,pval] = corr(D1',D2','Type','Spearman')
['end'];

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

