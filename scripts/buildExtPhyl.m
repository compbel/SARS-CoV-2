function [AM,badSeq,err,idsIntern] = buildExtPhyl(M,repSeq,repMuts,parMutPerSeq)
n = size(M,1);
m = size(M,2);
idsIntern = (cellstr(num2str((0:m)')))';
nRepSeq = length(repSeq);
mRepMut = length(repMuts);
Mext = zeros(n,m+mRepMut);
err = 0;
badSeq = [];
            Mext(1:n,1:m) = M;
            for j = 1:mRepMut
                ind = repSeq(find(parMutPerSeq == repMuts(j)));
                Mext(ind,m+j) = 1;
                Mext(ind,repMuts(j)) = 0;
                idsIntern{end+1} = [int2str(repMuts(j)) 'r'];
            end
            AM = buildPerfPhyl(Mext);
            G = digraph(AM);
            Mtree = zeros(n,m);
            for j = 1:n
                path = shortestpath(G,1,m+1+mRepMut + j);
                treeMut = sort(path(2:(end-1)));
                for k = treeMut
                    if k <= m+1
                        Mtree(j,k-1) = 1;
                    else
                        Mtree(j,repMuts(k-m-1)) = 1;
                    end
                end
                ind = find(Mtree(j,:) ~= M(j,:));
                if ~isempty(ind)
                    badSeq = [badSeq j];
                    err = err + length(ind);
                end
            end
            obj = 0;
            for u = 1:n
                for v = 1:m
                    if M(u,v) == 0
                        obj = obj + Mtree(u,v);
                    else
                        obj = obj - Mtree(u,v);
                    end
                end
            end
            ['aaa'];