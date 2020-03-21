function [AM,M,repMut] = extendTree(Mperf,Mnperf,Mobs,perfSeq,perfMut,nonPerfSeq,nonPerfMut,repMutnperf)
n = length(perfSeq) + length(nonPerfSeq);
m = length(perfMut) + length(nonPerfMut);
mRepMut = length(repMutnperf);
repMut = nonPerfMut(repMutnperf);
M = zeros(n,m+mRepMut);
M(perfSeq,1:m) = Mobs(perfSeq,:);
M(nonPerfSeq,[nonPerfMut (m+1):(m+mRepMut)]) = Mnperf;
M(nonPerfSeq,perfMut) = Mobs(nonPerfSeq,perfMut);

AM = buildPerfPhyl(M);
