function [mutGeneLabel,mutAALabel,nonsyn] = synMut(seq,ids,shift,varpos,ref)
geneSeg = [shift+1 265; 266 805; 806 2719; 2720 8554;  8555 10054; 10055 10972; 10973 11842; 11843 12091; 12092 12685; 12686 13024; 13025 13441; 13442 13468; 13468 16236; 16237 18039; 18040 19620; 19621 20658; 20659 21552; 13442 13480; 21563 25384; 25393 26220; 26245 26472; 26523 27191; 27202 27387; 27394 27759; 27756 27887; 27894 28259; 28274 29533; 29558 29674];
geneSeg = geneSeg - shift*ones(length(geneSeg),2);
geneName = {'5UTR','nsp1','nsp2','nsp3','nsp4','nsp5','nsp6','nsp7','nsp8','nsp9','nsp10','nsp12','nsp12','nsp13','nsp14','nsp15','nsp16','nsp11','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10'};
nGenes = length(geneSeg);
seqAA = cell(1,nGenes);
for i = 1:nGenes
    i
    gene(i).start = geneSeg(i,1);
    gene(i).end = geneSeg(i,2);
    gene(i).name = geneName{i};
    seqAA{i} = nt2aa(seq(:,gene(i).start:gene(i).end));
%     seqFile = ['testSyn' int2str(i) '.fas'];
%     if ~isfile(seqFile)
%         fastawrite(seqFile,ids,seq(:,gene(i).start:gene(i).end));
%     end
end

nonsyn = zeros(1,length(varpos));
geneID = zeros(1,length(varpos));
mutGeneLabel = cell(1,length(varpos));
mutAALabel = cell(1,length(varpos));
for i = 1:length(varpos)
    mut = varpos(i);
    mutAALabel{i} = '-';
    mutGeneLabel{i} = '-';
    g = intersect(find(geneSeg(:,1)<=mut),find(geneSeg(:,2)>=mut));
    if isempty(g)
        continue;
    end
    geneID(i) = g;
    mutGeneLabel{i} = gene(g).name; 
    if g == 1 
        continue;
    end
    mut1 = mut-gene(g).start+1;
    mutAA = ceil(mut1/3);
    un = unique(seqAA{g}(:,mutAA));    
    if length(un) > 1
        nonsyn(i) = 1;
        major = mode(seqAA{g}(:,mutAA));
        minor = setdiff(un,major);
        mutAALabel{i} = [major int2str(mutAA) minor];
    end
end

nMutGene = zeros(1,nGenes);
nNonsymGene = zeros(1,nGenes);
percNonsynGene = zeros(1,nGenes);
for i = 1:nGenes
    ind = find(geneID == i);
    nMutGene(i) = length(ind);
    nNonsymGene(i) = sum(nonsyn(ind),2);
    percNonsynGene(i) = sum(nonsyn(ind),2)/length(ind);
end
nonsyn = find(nonsyn);

% nposSyn = 0;
% nposNonSyn = 0;
% alph = 'ACTG';
% for i = 2:nGenes
%    i
%    refG = ref(gene(i).start:gene(i).end);
%    refGAA = nt2aa(refG);
%    for j = 1:size(refG,2)
%        jAA = ceil(j/3);
%        for k = 1:4
%            nuc = alph(k);
%            if nuc == refG(j)
%                continue;
%            end
%            varG = refG;
%            varG(j) = nuc;
%            varGAA = nt2aa(varG);
%            if varGAA(jAA) == refGAA(jAA)
%                nposSyn = nposSyn + 1;
%            else
%                nposNonSyn = nposNonSyn + 1;
%            end
%        end
%    end
% end
% nsyn = length(varpos)-length(nonsyn);
% p = nposSyn/(nposSyn+nposNonSyn)
% binopdf(nsyn,length(varpos),p)
% length(nonsyn)/(length(varpos)-length(nonsyn))
% ['end'];
% % seq = seq(:,gene(5).start:gene(5).end);
% % fastawrite('testSyn.fas',ids,seq);
