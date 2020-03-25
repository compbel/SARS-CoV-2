function AMNet = getTransNet(AMP,G,ids,times,seqMetadata,probHapAddr)
if ~isempty(probHapAddr)
    %imported from bootstrap
    load(probHapAddr);
end
n = length(ids);
m = size(AMP,1) - n - 1;
parLeafs = zeros(1,n);
for i = 1:n
    parLeafs(i) = find(AMP(:,m+1+i));
end

idsReg = ids;
for i = 1:n
    if strcmp(ids{i},'Zhejiang') || strcmp(ids{i},'Jiangsu') || strcmp(ids{i},'Hefei') || strcmp(ids{i},'Hangzhou') || strcmp(ids{i},'Anhui') || strcmp(ids{i},'Fujian')
        idsReg{i} = 'China';
    end
    if strcmp(ids{i},'Guangdong') || strcmp(ids{i},'Foshan') || strcmp(ids{i},'Jiangxi') || strcmp(ids{i},'Beijing') || strcmp(ids{i},'Sichuan') || strcmp(ids{i},'Yunnan')
        idsReg{i} = 'China';
    end
    if strcmp(ids{i},'Chongqing') || strcmp(ids{i},'Shenzhen') || strcmp(ids{i},'Shandong') || strcmp(ids{i},'Jingzhou') || strcmp(ids{i},'Tianmen') || strcmp(ids{i},'Guangzhou')
        idsReg{i} = 'China';
    end
    if strcmp(ids{i},'Nonthaburi')
        idsReg{i} = 'Thailand';
    end
    if strcmp(ids{i},'Czech-Republic')
        idsReg{i} = 'Czechia';
    end
    if strcmp(ids{i},'Sydney')
        idsReg{i} = 'Australia';
    end
    if strcmp(ids{i},'Korea')
        idsReg{i} = 'SouthKorea';
    end
end

haps = unique(parLeafs);
nHaps = length(haps);

H = transclosure(G);
AMcl = adjacency(H);
AMcl = AMcl(haps,haps);
H = transreduction(digraph(AMcl));
AMNet = adjacency(H);

hapsID = cell(1,nHaps);
hapsFirstTimes = [];
hapsLastTimes = [];
formatDate = 'dd/mm/yyyy';
countries = [];
for i = 1:nHaps
    ind = (parLeafs == haps(i));
    chHap = unique(idsReg(ind));
    countries = [countries chHap];
    infoHap = seqMetadata(ind);
%     infoHap = extractfield(seqMetadata(ind),"info");
    hapsID{i} = [chHap{1}];
    for j = 2:length(chHap)
        hapsID{i} = [hapsID{i} '/' chHap{j}];
    end
    hapStruct(i).Id = i;
    hapStruct(i).Sequences = find(ind);
    hapStruct(i).Label = hapsID{i};
    hapStruct(i).Count = length(ids(ind));
    hapStruct(i).Metadata = infoHap;
    t = min(times(parLeafs == haps(i)));
    tlast = max(times(parLeafs == haps(i)));
    hapsFirstTimes = [hapsFirstTimes t];
    hapsLastTimes = [hapsLastTimes tlast];
    if ~isnat(t)
        hapStruct(i).Date = datestr(t,formatDate);
        hapStruct(i).LastDate = datestr(tlast,formatDate);
    else
        hapStruct(i).Date = '';
        hapStruct(i).LastDate = '';
    end
end
% nViol = 0;
% timeDistNeigh = [];
% timeDistNeigh1 = [];
% for i = 1:nHaps
%     p = find(AMNet(:,i));
%     if isempty(p)
%         continue;
%     end
%      if ~isnat(hapsFirstTimes(i)) && ~isnat(hapsFirstTimes(p))
%          if (hapsFirstTimes(p) <= hapsFirstTimes(i)) && (hapsFirstTimes(i) <= hapsLastTimes(p))
%              timeDistNeigh = [timeDistNeigh 0];
%          else
%             timeDistNeigh = [timeDistNeigh days(hapsFirstTimes(i) - hapsLastTimes(p))];
%          end
%          tm1 = days(times(parLeafs == haps(p)) -  datetime(2019,12,1));
%          tm2 = days(times(parLeafs == haps(i)) -  datetime(2019,12,1));
%          timeDistNeigh1 = [timeDistNeigh1 min(min(abs(tm2 - tm1')))];
%      end
% %     if isnat(hapsTimes(i)) && ~isnat(hapsTimes(p))
% %         hapsTimes(i) = hapsTimes(p) + days(7);
% %         hapStruct(i).Date = datestr(hapsTimes(i),formatDate);
% %     end
%     if hapsFirstTimes(p) > hapsFirstTimes(i)
%         if days(hapsFirstTimes(p) - hapsFirstTimes(i)) >= 4
%             nViol = nViol + 1;
%         end
%         hapsFirstTimes(i) = hapsFirstTimes(p);
%         hapStruct(i).Date = datestr(hapsFirstTimes(i),formatDate);
%     end
% end
% 
% timeDistAll = [];
% for i = 1:nHaps
%     for p = (i+1):nHaps
%          if ~isnat(hapsFirstTimes(i)) && ~isnat(hapsFirstTimes(p))
%              tm1 = days(times(parLeafs == haps(p)) -  datetime(2019,12,1));
%              tm2 = days(times(parLeafs == haps(i)) -  datetime(2019,12,1));
%              timeDistAll = [timeDistAll min(min(abs(tm2 - tm1')))];
%          end
%     end
% end
% 
% [Nneigh,edges] = histcounts(timeDistNeigh1);
% [Nall,edges] = histcounts(timeDistAll);
% histogram(timeDistNeigh1)
% xlabel('time(days)')
% ylabel('number of potential transmission pairs')
% mean(timeDistNeigh1)
% confInterv(timeDistNeigh1)
% mean(timeDistAll)
% confInterv(timeDistAll)
% [h,p] = kstest2(timeDistNeigh1,timeDistAll)
% 
% 
% writetable(struct2table(hapStruct), 'transNetAttrN.csv');
% [s,t] = find(AMNet);
% for i = 1:length(s)
%     edge(i).Source = s(i);
%     edge(i).Target = t(i);
%     if ~isempty(probHapAddr)
%         edge(i).Weight = Probhap(s(i),t(i));
%     end
% end
% writetable(struct2table(edge), 'transNetAttrE.csv');

epid = zeros(1,nHaps);
for i = 1:nHaps
    md = hapStruct(i).Metadata;
    info = extractfield(md,'info');
    for j = 1:length(info)
        if contains(info{j},'Italy')
            epid(i) = 1;
        end
    end
end
epid = find(epid);

for i = 1:length(countries)
    if startsWith(countries{i},'USA')
        countries{i} = 'USA';
    end
    if startsWith(countries{i},'Netherlands')
        countries{i} = 'Netherlands';
    end
end
countries = unique(countries);
% for i = 1:length(countries)
%     H = digraph(AMNet);
%     figure
%     h = plot(H,'NodeLabel',hapsID,'Layout','force','NodeColor','b','MarkerSize',7);
%     h.NodeFontSize = 7;
%     ind = find(contains(hapsID,countries{i}));
%     highlight(h,ind,'NodeColor','r');
% end

H = digraph(AMNet);
figure
h = plot(H,'NodeLabel',hapsID,'Layout','force','NodeColor','b','MarkerSize',7);
h.NodeFontSize = 7;
ind = find(contains(hapsID,'Italy'));
highlight(h,ind,'NodeColor','r');
info = extractfield(hapStruct,'Metadata')';
info = info(ind);
% highlight(h,epid,'NodeColor','r');


% outdeg = sum(full(AMNet),2);
% outdd = histc(outdeg, unique(outdeg));
% k = find(outdd == 1,1);
% coeffs = polyfit(1:k,log(outdd(1:k))', 1)
% figure
% plot(1:k,log(outdd(1:k))');


Dist = distances(H);
D1 = [];
D2 = [];
hapsTimesNum = days(hapsFirstTimes - datetime(2019,12,1));
for i = 1:nHaps
    v = find(~isinf(Dist(i,:)));
    D1 = [D1 Dist(i,v)];
    D2 = [D2  (hapsTimesNum(v)-(hapsTimesNum(i)+2)*ones(1,length(v)))];
end
ind = find(~isnan(D2));
D1 = D1(ind);
D2 = D2(ind);
[rho,pval] = corr(D1',D2','Type','Spearman')
% scatter(D1,D2)
['end'];
