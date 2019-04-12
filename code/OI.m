function [Z1,k] = OI(A,iteration)
n = size(A);
n = n(1);
Z0 = eye(n);
k = 1;
while k<=iteration
    Y = A*Z0;
    [Z1,R] = qr(Y);
    k = k+1;
    if norm(Z1-Z0) < 1e-8
        break
    end
    Z0 = Z1;
end
k = k-1;
