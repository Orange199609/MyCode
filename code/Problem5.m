n = input('��������ά��:');
lambda = sort(10*rand(n,1));  %�������һ��Ԫ�ش�С��0-10֮���nά������ΪA������ֵ
x = randn(n);
A = x*diag(lambda)*inv(x);   %���ɴ��ֽ����
iter_limit = 1000;
fprintf('Eigenvalue randomly generated are:\n');
lambda
[z1,k1] = OI(A,iter_limit);
A1 = z1'*A*z1;   %��z1�ԶԽǾ�������ԭ
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