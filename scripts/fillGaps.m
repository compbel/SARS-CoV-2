function seqNew = fillGaps(seq)
nucs = 'ACTG';
seqNew = seq;
for i = 1:size(seq,1)
    for j = 1:size(seq,2)
        if ~contains(nucs,seq(i,j))
            freq = zeros(1,4);
            for k = 1:4
                freq(k) = sum(seq(:,j) == nucs(k));
            end
            l = find(freq == max(freq));
            l = l(1);
            seqNew(i,j) = nucs(l);
        end
    end
end