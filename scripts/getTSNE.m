function [] = getTSNE(seq,ids,varpos)

n = length(ids);
l = size(seq,2);
idsReg = ids;
for i = 1:n
    if strcmp(ids{i},'Zhejiang') || strcmp(ids{i},'Jiangsu') || strcmp(ids{i},'Hefei') || strcmp(ids{i},'Hangzhou') || strcmp(ids{i},'Anhui') || strcmp(ids{i},'Fujian')
        idsReg{i} = 'China';
    end
    if strcmp(ids{i},'Guangdong') || strcmp(ids{i},'Foshan') || strcmp(ids{i},'Jiangxi') || strcmp(ids{i},'Beijing') || strcmp(ids{i},'Sichuan') || strcmp(ids{i},'Yunnan')
        idsReg{i} = 'China';
    end
    if strcmp(ids{i},'Chongqing') || strcmp(ids{i},'Shenzhen') || strcmp(ids{i},'Shandong') || strcmp(ids{i},'Jingzhou') || strcmp(ids{i},'Tianmen')
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
        idsReg{i} = 'South-Korea';
    end
end

[haps,unI,hapI] = unique(seq,'rows');
nHaps = size(haps,1);
hapsID = cell(1,nHaps);
for i = 1:nHaps
    ind = find(hapI == i);
    names = unique(idsReg(ind));
    hapsID{i} = [names{1}];
    for j = 2:length(names)
        hapsID{i} = [hapsID{i} '/' names{j}];
    end
end


X = zeros(nHaps,l);
for i = 1:nHaps
    for j = 1:l
        if haps(i,j) == 'A'
            X(i,j) = 1;
        end
        if haps(i,j) == 'C'
            X(i,j) = 2;
        end
        if haps(i,j) == 'T'
            X(i,j) = 3;
        end
        if haps(i,j) == 'G'
            X(i,j) = 4;
        end
    end
end

% optVal = 0;
% optY = [];
% optEva = [];
% for i = 1:40
%     i
%     Y = tsne(X,'Algorithm','exact','Distance','hamming');
%     eva = evalclusters(Y,'kmeans','gap','KList',[1:8]);
%     if max(eva.CriterionValues) > optVal
%         optVal = max(eva.CriterionValues);
%         optY = Y;
%         optEva = eva;
%     end
% end
% figure
% plot(optEva)
% k = optEva.OptimalK;

load('clust2.mat');
k = 5;
% export_fig 'clust2eva.svg' -m5

idx = kmeans(optY,k,'MaxIter',1000,'Replicates',100);
% idx = spectralcluster(Y,k);
figure
gscatter(optY(:,1),optY(:,2),idx,'bgmrk',[],30,'doleg','off')
% labelpoints (optY(:,1)+0.1,optY(:,2),hapsID,'N', 0.5,'FontSize',7)
text(optY(:,1)+0.1,optY(:,2),hapsID,'FontSize',7)
% figure
% textscatter(optY(:,1),optY(:,2),hapsID,'TextDensityPercentage',100)

% strOut = [num2cell(optY(:,1)) num2cell(optY(:,2)) num2cell(idx) hapsID'];
% writetable(cell2table(strOut),'tsne.csv');
% 
% haps = haps(:,varpos);
k = 4;
% kmerAll = nmercount(haps(1,:), k);
% kmerAll = kmerAll(:,1);
% for i = 2:nHaps
%     kmer = nmercount(haps(i,:), k);
%     kmerAll = union(kmerAll,kmer(:,1));
% end
% nkm = length(kmerAll);
% X = zeros(nHaps,nkm);
% for i = 1:nHaps
%    kmer = nmercount(haps(i,:), k);
%    for j = 1:length(kmer)
%        ind = find(strcmp(kmerAll,kmer{j,1}));
%        X(i,ind) = kmer{j,2};
%    end
% end
% Y = tsne(X,'Algorithm','exact','Distance','correlation');
% k = 5;
% idx = kmeans(Y,k,'MaxIter',1000,'Replicates',100);
% % idx = spectralcluster(Y,k);
% figure
% gscatter(Y(:,1),Y(:,2),idx,'bgmrk',[],30)
% text(Y(:,1)+0.1,Y(:,2),hapsID,'FontSize',7)

['end'];