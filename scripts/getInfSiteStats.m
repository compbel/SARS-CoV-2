function [infSiteStats,genotMatr,minorMinor] = getInfSiteStats(hapAnalize,freq)

        n = size(hapAnalize,1);
        m = size(hapAnalize,2);
        consSeq = '';
        for i = 1:m
            major = '';
            currMaxFreq = 0;
            nucs = unique(hapAnalize(:,i));
            freqAlleles = zeros(1,length(nucs));
            for j = 1:length(nucs)
                ind = (hapAnalize(:,i) == nucs(j));
                freqAlleles(j) = sum(freq(ind),2);
                if freqAlleles(j) > currMaxFreq
                    currMaxFreq = freqAlleles(j);
                    major = nucs(j);
                end
            end
            consSeq = [consSeq major];
        end
        
        genotMatr = zeros(n,m);
        for i = 1:n
            for j = 1:m
                if hapAnalize(i,j) ~= consSeq(j)
                    genotMatr(i,j) = 1;
                end
            end
        end
        
        infSiteStats = 0;
        minorMinor = [];
        for i = 1:m
            for j = (i+1):m
                if (sum(genotMatr(:,i) + genotMatr(:,j) == 0,1) > 0) && (sum(genotMatr(:,i) - genotMatr(:,j) > 0,1) > 0) && (sum(genotMatr(:,j) - genotMatr(:,i) > 0,1) > 0) && (sum(genotMatr(:,i).*genotMatr(:,j) > 0,1) > 0)
                    infSiteStats = infSiteStats + 1;
                    ind = find(genotMatr(:,i).*genotMatr(:,j) > 0);
                    minorMinor = [minorMinor; [ind i*ones(length(ind),1) j*ones(length(ind),1)]];
                end
            end
        end
        
        
        ['End']