function [leftDesc,rightDesc] = getDesc(bintree)
% bintree: [leftChild righChild parent label frequency timet fitness]
nVert = size(bintree,1);
leftDesc = cell(1,nVert);
rightDesc = cell(1,nVert);

AM = treeAMdirect(bintree);
G = digraph(AM);
sort = flip(toposort(G));

for v = sort
    if bintree(v,1) + bintree(v,2) ~= 0
        leftDesc{v} = [leftDesc{bintree(v,1)} rightDesc{bintree(v,1)} bintree(v,1)];
        rightDesc{v} = [leftDesc{bintree(v,2)} rightDesc{bintree(v,2)} bintree(v,2)];
    end
end
