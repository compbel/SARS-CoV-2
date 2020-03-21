function [seqAA,allMutNames] = getAA(seq,shift)
geneSeg = [shift+1 265; 266 805; 806 2719; 2720 8554;  8555 10054; 10055 10972; 10973 11842; 11843 12091; 12092 12685; 12686 13024; 13025 13441; 13442 13468; 13468 16236; 16237 18039; 18040 19620; 19621 20658; 20659 21552; 21563 25384; 25393 26220; 26245 26472; 26523 27191; 27202 27387; 27394 27759; 27756 27887; 27894 28259; 28274 29533; 29558 29674];
geneSeg = geneSeg - shift*ones(length(geneSeg),2);
geneName = {'5UTR','nsp1','nsp2','nsp3','nsp4','nsp5','nsp6','nsp7','nsp8','nsp9','nsp10','nsp12','nsp12','nsp13','nsp14','nsp15','nsp16','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10'};
nGenes = length(geneSeg);
seqAA = [];
allMutNames = {};
for i = 1:nGenes
    i
    seqAAg = nt2aa(seq(:,geneSeg(i,1):geneSeg(i,2)));
    seqAA = [seqAA seqAAg];
    nm = cell(1,size(seqAAg,2));
    nm(:) = {geneName{i}};
    allMutNames = [allMutNames nm];
    
%     n = size(seqAA,1);
%     seqFile = ['testSyn_' geneName{i} '.fas'];
%     if ~isfile(seqFile)
%         fastawrite(seqFile,(cellstr(num2str((1:n)')))',seqAAg);
%     end
end
