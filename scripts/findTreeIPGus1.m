function [AM,M_infer_full,err,objval] = findTreeIPGus1(M,repMut,cutoff)
eps = 0.01;

% repMut = [repMut repMut];
% repmut = 1:size(M,2);

M_full = M;
[M, mut_red, mut_full] = unique(M_full', 'rows');
M = M';
repMut_full = repMut;
repMut = mut_full(repMut)';
[M, seq_red, seq_full] = unique(M, 'rows');

freqMut = zeros(1,length(mut_red));
for i = 1:length(mut_red)
    freqMut(i) = sum(mut_full == i,1);
end

freqSeq = zeros(1,length(seq_red));
for i = 1:length(seq_red)
    freqSeq(i) = sum(seq_full == i,1);
end

n = size(M,1);
mOrig = size(M,2);
mRepMut = length(repMut);
m = mOrig+mRepMut;
allMuts = [1:mOrig repMut];

fileID = fopen('ILP.lp','w');
fprintf(fileID,'%s\n','minimize');
for i = 1:n
    for j = 1:mOrig
        ind = find(allMuts == j);
        coeff = [int2str(freqSeq(i)*freqMut(j)) ' '];
        if M(i,j) == 0
            for j1 = ind
                Yij1 = ['Y' int2str(i) ',' int2str(j1)];
                fprintf(fileID,'%s',['+ ' coeff Yij1 ' ']);
            end
        else
            for j1 = ind
                Yij1 = ['Y' int2str(i) ',' int2str(j1)];
                fprintf(fileID,'%s',['- ' coeff Yij1 ' ']);
            end
        end
    end
end

fprintf(fileID,'\n%s\n',' st');

for p = 1:m
    for q = (p+1):m
        [p q]
         D00 = ['D' int2str(p) ',' int2str(q) ',0,0'];
         D01 = ['D' int2str(p) ',' int2str(q) ',0,1'];
         D10 = ['D' int2str(p) ',' int2str(q) ',1,0'];
         D11 = ['D' int2str(p) ',' int2str(q) ',1,1'];
        for i = 1:n
            Yip = ['Y' int2str(i) ',' int2str(p)];
            Yiq = ['Y' int2str(i) ',' int2str(q)];
            cons = [D11 ' - ' Yip ' - ' Yiq ' >= -1'];
            fprintf(fileID ,'%s\n',cons);
            cons = [D00 ' + ' Yip ' + ' Yiq ' >= 1'];
            fprintf(fileID ,'%s\n',cons);
            cons = [D10 ' - ' Yip ' + ' Yiq ' >= 0'];
            fprintf(fileID ,'%s\n',cons);
            cons = [D01 ' + ' Yip ' - ' Yiq ' >= 0'];
            fprintf(fileID ,'%s\n',cons);
        end
        cons = [D00 ' + ' D01 ' + ' D10 ' + ' D11 ' <= 3'];
        fprintf(fileID ,'%s\n',cons);
    end
end

for j = 1:m
    ind = find(allMuts == j);
    if length(ind) >= 2
        for i = 1:n
            cons = [];
            for p = ind
                Yip = ['Y' int2str(i) ',' int2str(p)];
                cons = [cons '+ ' Yip ' '];
            end
            cons = [cons ' <= 1'];
            fprintf(fileID ,'%s\n',cons);
        end
    end
end

fprintf(fileID,'%s\n','binaries');
for p = 1:m
    for q = (p+1):m
         [p q]
         D00 = ['D' int2str(p) ',' int2str(q) ',0,0'];
         D01 = ['D' int2str(p) ',' int2str(q) ',0,1'];
         D10 = ['D' int2str(p) ',' int2str(q) ',1,0'];
         D11 = ['D' int2str(p) ',' int2str(q) ',1,1'];
         fprintf(fileID,'%s\n',D00);
         fprintf(fileID,'%s\n',D01);
         fprintf(fileID,'%s\n',D10);
         fprintf(fileID,'%s\n',D11);
    end
end
for i = 1:n
    for j = 1:m
        [n m]
        Yij = ['Y' int2str(i) ',' int2str(j)];
        fprintf(fileID,'%s\n',Yij);
    end
end
fprintf(fileID,'%s\n','end');
fclose(fileID);

% Mgus = 2*ones(n,m);
% writematrix(Mgus);
% [results,status] = perl('bbmisstest.pl','Mgus.txt','ILPGus.lp');
% system('copy objective+constraints ILP.lp');
clear model;
model = gurobi_read('ILP.lp');

params.outputflag = 1;
params.timelimit = 2500;
% params.NodefileStart = 30;
% params.MIPFocus = 2;
params.threads = 18;
params.mipgap = 0.0075;
% params.cutoff = cutoff;
result = gurobi(model, params);
if strcmp(result.status,'CUTOFF')
    AM = [];
    M_infer_full = [];
    err = Inf;
    objval = Inf;
    return;
end
objval = result.objval;

matNames = model.varnames(1:(n*m));
M_infer = zeros(n,m);
for i = 1:n
    for j = 1:m
        Yij = ['Y' int2str(i) ',' int2str(j)];
        ind = find(ismember(matNames,Yij));
        M_infer(i,j) = result.x(ind);
    end
end


%     m = size(M,2);
%     AM = buildPerfPhyl(M_infer);
%     G = digraph(AM);
% %     idsIntern = (cellstr(num2str((0:m)')))';
% %     for j = 1:mRepMut
% %         idsIntern{end+1} = [int2str(repMut(j)) 'r'];
% %     end
% %     idsLeafs = (cellstr(num2str((1:n)')))';
%     idsIntern = (cellstr(num2str([0 mut_red']')))';
%     for j = 1:mRepMut
%         idsIntern{end+1} = [int2str(repMut_full(j)) 'r'];
%     end
%     idsLeafs = (cellstr(num2str((1:n)')))';
%     mutMarkers = cell(1,m+1+mRepMut);
%     mutMarkers(:) = {'o'};
%     patMarkers = cell(1,n);
%     patMarkers(:) = {'s'};
%     figure
%     h = plot(G,'NodeLabel',[idsIntern idsLeafs],'Layout','layered','Marker',[mutMarkers patMarkers],'NodeColor','b','MarkerSize',7);
%     highlight(h,(m+mRepMut+2):(m+1+mRepMut+n),'NodeColor','r');
%     highlight(h,(m+2):(m+1+mRepMut),'NodeColor','g');
    
M_infer_full = zeros(size(M_full,1),size(M_full,2)+mRepMut);
M_infer = M_infer(seq_full,:);
for j = 1:size(M_full,2)
    M_infer_full(:,j) = M_infer(:,mut_full(j));
end
for j = 1:mRepMut
    M_infer_full(:,size(M_full,2) + j) = M_infer(:,size(M,2)+j);
end

m = size(M_full,2);
n = size(M_full,1);
AM = buildPerfPhyl(M_infer_full);
G = digraph(AM);
Mtree = zeros(n,m);
err = 0;
for j = 1:n
    path = shortestpath(G,1,m+1+mRepMut + j);
    treeMut = sort(path(2:(end-1)));
    for k = treeMut
        if k <= m+1
            Mtree(j,k-1) = 1;
        else
            Mtree(j,repMut_full(k-m-1)) = 1;
        end
    end
    ind = find(Mtree(j,:) ~= M_full(j,:));
    if ~isempty(ind)
        err = err + length(ind);
    end
end
    

    
['end']
end
    

