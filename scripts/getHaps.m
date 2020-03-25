function hapsID = getHaps(AM,G,ids,times,seqMetadata)
n = length(ids);
m = size(AM,1) - n - 1;
parLeafs = zeros(1,n);
for i = 1:n
    parLeafs(i) = find(AM(:,m+1+i));
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
hapsTimes = cell(1,nHaps);
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
    hapStruct(i).Label = hapsID{i};
    hapStruct(i).Count = length(ids(ind));
    hapStruct(i).Metadata = infoHap;
    t = min(times(parLeafs == haps(i)));
    hapsTimes{i} = t;
    if ~isnat(t)
        hapStruct(i).Date = datestr(t,formatDate);
    else
        hapStruct(i).Date = '';
    end
end

