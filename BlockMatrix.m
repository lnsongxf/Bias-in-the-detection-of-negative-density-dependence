function    B = BlockMatrix(A,dx)

[M,N]=size(A);

n=M/dx;
m=N/dx;
B = zeros(n,m);
for i=1:n
    for j=1:m
        temp = A((i-1)*dx+1:i*dx,(j-1)*dx+1:j*dx);
        B(i,j)=sum(temp(:));
    end
end


