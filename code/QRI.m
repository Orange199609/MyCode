function [A1,k] = QRI(A,iteration)
k = 1;
A0 = A;
while k <= iteration
    [Q,R] = qr(A0);
    A1 = R*Q;
    k = k+1;
    if norm(A1-A0) < 1e-8
        break
    end
    A0 = A1;
end
k = k-1;