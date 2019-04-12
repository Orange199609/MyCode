function [error1,error2] = udsys(A,b)
% Method 1
matrix1 = inv(A'*A);
matrix2 = A'*b;
x1 = matrix1*matrix2;
error1 = norm(A*x1-b);
x2 = A\b;
error2 = norm(A*x2-b);
fprintf('Error of the 1st method (using A_transpose*A) is %.4f\n',error1);
fprintf('Error of the 2nd method (using operator) is %.4f\n',error2);
