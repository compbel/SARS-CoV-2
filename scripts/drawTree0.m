function [] = drawTree(AM,n,m,repMuts,patID,varpos,shift,mutLabel,nonsyn)

    varpos = varpos + shift;
    G = digraph(AM);
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
    nonsyn = nonsyn + ones(1,length(nonsyn));
    figure
    h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
    highlight(h,(m+mRepMut+2):(m+1+mRepMut+n),'NodeColor','r');
%     highlight(h,nonsym,'NodeColor','y');
%     highlight(h,[repMuts + 1 (m+2):(m+1+mRepMut)],'NodeColor','g');
    
    mutFreq = sum(AM(1:(m+mRepMut+1),(m+mRepMut+2):end),2)+1;
    unobs = find(mutFreq == 1);
    mutFreq = mutFreq/(sum(mutFreq,1));
    markSizes = 700*mutFreq;
    Gmut = digraph(AM(1:(m+mRepMut+1),1:(m+mRepMut+1)));
    figure
%     h = plot(Gmut,'Layout','layered','NodeColor','b','MarkerSize',markSizes,'LineWidth', 2,'NodeLabel',idsIntern);
    h = plot(Gmut,'Layout','layered','NodeColor','b','LineWidth', 2,'NodeLabel',idsIntern,'MarkerSize',7);
    highlight(h,nonsyn,'NodeColor','r');
%     highlight(h,unobs,'NodeColor','w');
%     highlight(h,[repMuts + 1 (m+2):(m+1+mRepMut)],'NodeColor','g');
['end'];