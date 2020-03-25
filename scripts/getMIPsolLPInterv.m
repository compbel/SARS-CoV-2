function [M_infer,x,u,v,t] = getMIPsolLPInterv(var,n,m)

varInd = [(n+m)*(n+m) n+m n+m (n+m)*(n+m-1)/2];

ivar = 1;

x = zeros(n+m,n+m);
for i = 1:(n+m)
    for j = 1:(n+m)
%         if i == j
%             continue;
%         end
        x(i,j) = var(ivar);
        ivar = ivar + 1;
    end
end

u = zeros(1,n+m);
v = zeros(1,n+m);
for i = 1:(n+m)
    u(i) = var(ivar);
    ivar = ivar + 1;
end
for i = 1:(n+m)
    v(i) = var(ivar);
    ivar = ivar + 1;
end

t = zeros(n+m,n+m);
for i=1:(n+m)
    for j = (i+1):(n+m)
        if var(ivar) == 1
            t(i,j) = 1;
        else
            t(j,i) = 1;
        end
        ivar = ivar + 1;
    end
end

M_infer = x((n+1):(n+m),1:n)';


