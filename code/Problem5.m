n = input('创建方阵维度:');
lambda = sort(10*rand(n,1));  %随机生成一个元素大小在0-10之间的n维向量作为A的特征值
x = randn(n);
A = x*diag(lambda)*inv(x);   %生成待分解矩阵
iter_limit = 1000;
fprintf('Eigenvalue randomly generated are:\n');
lambda
[z1,k1] = OI(A,iter_limit);
A1 = z1'*A*z1;   %用z1对对角矩阵做还原
error1 = norm(lambda-sort(diag(A1)));
fprintf('Eigenvalue calculate by OI are:\n');
sort(diag(A1))
fprintf('The error of OI result is %.4f\n',error1);
fprintf('OI process took %.1f steps\n',k1);
[z2,k2] = QRI(A,iter_limit);
error2 = norm(lambda-sort(diag(z2)));
fprintf('Eigenvalue calculated by QRI are:\n');
sort(diag(z2))
fprintf('The error of QRI result is %.4f\n',error2);
fprintf('QRI process took %.1f steps\n',k2);