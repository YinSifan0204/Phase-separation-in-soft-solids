function D1 = chebDiffmat(N,L)
%Derivative matrix of first order
%  

D1  =zeros(N+1,N+1);
    for i = 1:N+1
        for j = i+1 : 2 : N+1
              D1(i,j) = 2*(j-1);
        end
    end
    D1(1,:) = D1(1,:)/2;
    D1 = 2*D1/L;