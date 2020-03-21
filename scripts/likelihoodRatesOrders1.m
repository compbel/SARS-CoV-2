function f = likelihoodRatesOrders1(t,theta,mutOrders)
f = 0;
m = length(t)-1;
for i = 1:m
    for j = 2:length(mutOrders{i})
        u = mutOrders{i}(j-1);
        v = mutOrders{i}(j);
        timeInt = t(v) - t(u);
        f = f - log(theta(i)*timeInt) + theta(i)*timeInt;
    end
    f = f + theta(i)*(t(m+1) - t(mutOrders{i}(end)));
end


% if nargout > 1 % gradient required
%     gradf = zeros(1,m);
%     for i = 1:m
%         gradf(i) = theta(i) - deg(i)/t(i);
%     end
% end
