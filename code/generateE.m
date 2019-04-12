function [C,E] = generateE(n)
C = zeros(n,n);
for row = 1:n
    C(row,row) = row*(n-row+1);
    for col = row+1:n
        C(row,col) = C(row,col-1)-row;
    end
    for col = 1:row-1
        C(row,col) = C(col,row);
    end
end
E = (1/(n+1))*C;

