function [L,P] = myspchol(A)
n = size(A,1);
base = zeros(n);
iter_max = factorial(n);     %最多随机试验n的阶乘次（有Bootstrapping的思想）
sparse_min = n*n;            %记录目前最小的稀疏值
L = zeros(n);                %记录稀疏值最小情况下的分解结果
P = zeros(n);                %记录pivoting规则
for i = 1:iter_max
    col_index = randperm(n);
    P0 = base;
    for j = 1:n
        P0(j,col_index(j)) = 1;
    end
    L0 = chol(P0*A*P0','lower');
    if nnz(L0) < sparse_min
        L = L0;
        P = P0;
    end
end
end
