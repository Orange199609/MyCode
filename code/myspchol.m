function [L,P] = myspchol(A)
n = size(A,1);
base = zeros(n);
iter_max = factorial(n);     %����������n�Ľ׳˴Σ���Bootstrapping��˼�룩
sparse_min = n*n;            %��¼Ŀǰ��С��ϡ��ֵ
L = zeros(n);                %��¼ϡ��ֵ��С����µķֽ���
P = zeros(n);                %��¼pivoting����
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
