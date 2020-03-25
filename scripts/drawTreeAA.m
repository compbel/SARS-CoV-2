function [] = drawTreeAA(AM,n,m,repMuts,patID,varpos,allMutNamesAA)

    G = digraph(AM);
    mRepMut = length(repMuts);
    idsIntern = cellstr([{'0'} allMutNamesAA(varpos)]);
    for j = 1:mRepMut
        idsIntern(end+1) = allMutNamesAA(varpos(repMuts(j)));
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
    figure
    h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
    highlight(h,(m+mRepMut+2):(m+1+mRepMut+n),'NodeColor','r');

    
    mutFreq = sum(AM(1:(m+mRepMut+1),(m+mRepMut+2):end),2)+1;
    unobs = find(mutFreq == 1);
    mutFreq = mutFreq/(sum(mutFreq,1));
    markSizes = 700*mutFreq;
    Gmut = digraph(AM(1:(m+mRepMut+1),1:(m+mRepMut+1)));
    figure
    h = plot(Gmut,'Layout','layered','NodeColor','b','MarkerSize',markSizes,'LineWidth', 2,'NodeLabel',idsIntern);
%     h = plot(Gmut,'Layout','force','NodeColor','b','LineWidth', 2,'NodeLabel',idsIntern,'MarkerSize',7);

%     highlight(h,unobs,'NodeColor','w');
%     highlight(h,[repMuts + 1 (m+2):(m+1+mRepMut)],'NodeColor','g');
['end'];