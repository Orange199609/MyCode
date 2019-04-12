%Problem 2
function [L] = mychol(A)
%using Cholesky Decomposition to decomposite a Positive Definite Matrix
%A:Positive Definite Matrix
[m,n] = size(A);
L = zeros(m,n);
L(1,1) = sqrt(A(1,1));
for row=2:m
    L(row,1) = A(row,1)/L(1,1);
end

for k = 2:n
    L(k,k) = A(k,k);
    for col = 1:k-1
        L(k,k) = L(k,k) - L(k,col)^2;
    end
    L(k,k) = sqrt(L(k,k));
    for i = k+1:m
        L(i,k) = A(i,k);
        for j = 1:k-1
            L(i,k) = L(i,k) - L(i,j)*L(k,j);
        end
        L(i,k) = L(i,k)/L(k,k);
    end
end
