function [] = drawTree(AM,n,m,repMuts,patID,varpos,shift,mutLabel,nonsyn)

    varpos = varpos + shift;
    mRepMut = length(repMuts);
    if isempty(mutLabel)
        idsIntern = (cellstr(num2str([0 varpos]')))';
        for j = 1:mRepMut
            idsIntern{end+1} = [int2str(varpos(repMuts(j))) 'r'];
        end
    else
        idsIntern = [{'0'} mutLabel];
        for j = 1:mRepMut
            idsIntern{end+1} = [mutLabel{repMuts(j)}];
        end
    end
    if isempty(patID)
        idsLeafs = (cellstr(num2str((1:n)')))';
    else
        idsLeafs = patID;
    end
    mutMarkers = cell(1,m+1+mRepMut);
    mutMarkers(:) = {'o'};
    patMarkers = cell(1,n);
    patMarkers(:) = {'s'};
    if ~isempty(nonsyn)
        nonsyn = nonsyn + ones(1,length(nonsyn));
    end
    
    idsAll = [idsIntern idsLeafs];
    markersAll = [mutMarkers patMarkers];
    nonisol = union(find(sum(AM,2) > 0),find(sum(AM,1) > 0)); 
    haps = (length(nonisol)-n+1):(length(nonisol));
    AMnonisol = AM(nonisol,nonisol);
    idsAll = idsAll(nonisol);
    markersAll = markersAll(nonisol);
    G = digraph(AMnonisol);
    
    figure
    h = plot(G,'NodeLabel',idsAll,'Layout','layered','Marker',markersAll,'NodeColor','b','MarkerSize',7);
    highlight(h,haps,'NodeColor','r');
%     highlight(h,nonsym,'NodeColor','y');
%     highlight(h,[repMuts + 1 (m+2):(m+1+mRepMut)],'NodeColor','g');
    
    mutFreq = sum(AM(1:(m+mRepMut+1),(m+mRepMut+2):end),2)+1;
    unobs = find(mutFreq == 1);
    mutFreq = mutFreq/(sum(mutFreq,1));
    markSizes = 700*mutFreq;
    AMmut = AM(1:(m+mRepMut+1),1:(m+mRepMut+1));
    nonisol = union(find(sum(AMmut,2) > 0),find(sum(AMmut,1) > 0)); 
    AMmut = AMmut(nonisol,nonisol);
    idsIntern = idsIntern(nonisol);
    nonSynVect = logical(zeros(1,m+mRepMut+1));
    nonSynVect(nonsyn) = 1;
    nonSynVect = nonSynVect(nonisol); 
    Gmut = digraph(AMmut);
    figure
%     h = plot(Gmut,'Layout','layered','NodeColor','b','MarkerSize',markSizes,'LineWidth', 2,'NodeLabel',idsIntern);
    h = plot(Gmut,'Layout','layered','NodeColor','b','LineWidth', 2,'NodeLabel',idsIntern,'MarkerSize',7);
    if ~isempty(nonsyn)
        highlight(h,nonSynVect,'NodeColor','r');
    end
%     export_fig 'mutTree.png' -m5
%     highlight(h,unobs,'NodeColor','w');
%     highlight(h,[repMuts + 1 (m+2):(m+1+mRepMut)],'NodeColor','g');
['end'];