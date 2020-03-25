function times = nCoV_param(AMP,n,m,timeLeafs,l)

eps = 0.001;
eps1 = 0.001;
maxIter = 10;
expTheta = l*3*10^(-4);
minTheta = l*10^(-4);
maxTheta = l*6*10^(-4);

outdeg = sum(AMP,2);
muts = find(outdeg > 0);
people = find(outdeg == 0);

for i = 1:(m+1)
    stree(i).children = intersect(find(AMP(i,:) > 0),muts)';
    stree(i).parent = find(AMP(:,i) > 0);
    stree(i).time = 0;
    stree(i).rate = 0;
    stree(i).people = intersect(find(AMP(i,:) > 0),people)'-m-1;
end

dt1 = datetime([2020 1 1]);
dt2 = datetime([2020 12 31]);
timeLeafsRel = (timeLeafs - min(timeLeafs))/(dt2 - dt1);

internTimes = -ones(1,m+1);
for i = 1:(m+1)
    if ~isempty(stree(i).people)
        t = min(timeLeafsRel(stree(i).people));
        if ~isempty(t)
            internTimes(i) = t;
        end
    end
end

rates = expTheta*ones(1,m+1);
mutOrders = cell(1,m+1);
for i = 1:(m+1)
    mutOrders{i} = [i stree(i).children];
end
[times,likel] = findTimeFixedOrderRate5(stree,AMP(1:(m+1),1:(m+1)),mutOrders,rates,internTimes,minTheta,maxTheta);

['End']