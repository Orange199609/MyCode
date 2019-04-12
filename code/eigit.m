function [lam,u,iter,v] = eigit(A,tol)
% Solves EVP to determine dominant eigenvalue and associated vector
% Sample call: [lam u iter] = eigit(A,tol)
% A is a square matrix, tol is the accuracy
% lam is the dominant eigenvalue, u is the associated vector
% iter is the number of iterations required
n = size(A);
n = n(1);
err = 100*tol;
u0 = ones(n,1); iter = 0;
while err>tol
v = A*u0;
u1 = (1/max(v))*v;
err = max(abs(u1-u0));
u0 = u1; iter = iter+1;
end
u = u0; 
lam = max(v);
%lam = v(find(v1==min(abs(v1))));