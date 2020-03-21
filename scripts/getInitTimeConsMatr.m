function [currTime,EMOrders,A,rhs] = getInitTimeConsMatr(mutOrders,stree,minRate,maxRate,maxTime,m,leafs)
    nLeafs = length(leafs);
    theta0 = 0.5*(maxRate + minRate);
    currTime = zeros(1,m+1);
    AMOrders = zeros(m,m);
    EMOrders = zeros(m-1,2*m);
%     EMOrders = zeros(m-1+nLeafs,m+1);
    iEdge = 1;
        for i = 1:m
            for j = 2:length(mutOrders{i})
                AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
                EMOrders(iEdge,mutOrders{i}(j)) = 1;
                EMOrders(iEdge,mutOrders{i}(j-1)) = -1;
                iEdge = iEdge + 1;
            end
%             if ismember(i,leafs)
%                 EMOrders(iEdge,m+1) = 1;
%                 EMOrders(iEdge,mutOrders{i}(end)) = -1;
%                 iEdge = iEdge + 1;
%             end
        end
    G = digraph(AMOrders);
    % plot(G);
    % G = digraph(AM);
    v = bfsearch(G,1);
    for i = 1:m
        currTime(v(i)) = (i-1)/theta0;
    end

    currTime(m+1) = max(currTime(leafs));
    
    A = zeros(m-1+nLeafs,m+1);
    rhs = zeros(m-1+nLeafs,1);
    iEq = 0;
    for i = 1:m
        if ~isempty(stree(i).children)
            for j = 2:length(mutOrders{i})
                iEq = iEq + 1;
                A(iEq,mutOrders{i}(j-1)) = 1;
                A(iEq,mutOrders{i}(j)) = -1;
            end
        end
    end
    for i = leafs
        iEq = iEq + 1;
        A(iEq,i) = 1;
        A(iEq,m+1) = -1;
    end
    
    A(end+1,m+1) = 1;
    rhs(end+1) = maxTime;
    currTime = maxTime*(currTime/sum(currTime,2));
    % Abd = zeros(nIntern,m+1);
    % rhsbd = zeros(nIntern,1);
    % iEq = 0;
    % for i = 13
    %     iEq = iEq + 1;
    %     Abd(iEq,i) = 1;
    %     Abd(iEq,m+1) = -1;
    %     rhsbd(iEq) = -deg(i)/maxRate;
    %     iEq = iEq + 1;
    %     Abd(iEq,i) = -1;
    %     Abd(iEq,m+1) = 1;
    %     rhsbd(iEq) = deg(i)/minRate;
    % end
    % 
    % A = [A; Abd];
    % rhs = [rhs; rhsbd];