%Problem4 Script
for k = 3:6
    fprintf('µ±n=%.0f ±\n',k);
    A = generateA(k);
    A_SPD = A'*A;
    A_inv = inv(A);
    P = inv(A_SPD);
    Q = A_inv*A_inv';
    R = invhilb(k)*invhilb(k)';
    error1 = norm(P-R);
    error2 = norm(Q-R);
    fprintf('P error is %.4f\n',error1);
    fprintf('Q error is %.4f\n',error2);
end