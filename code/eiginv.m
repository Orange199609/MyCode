function [lam,u,iter] = eiginv(A,mu,tol)
% Determines eigenvalue of A closest to mu with a tolerance tol.
% Sample call: [lam u] = eiginv(A,mu,tol)
% lam is the eigenvalue and u the corresponding eigenvector.
n = size(A);
n = n(1);
err = 100*tol;
B = A-mu*eye(n,n);
u0 = ones(n,1);
iter = 0;
while err>tol
v = B\u0; f = 1/max(v);
u1 = f*v;
err = max(abs(u1-u0));
u0 = u1; iter = iter+1;
end
u = u0; lam = mu+f;