function D2 = chebDiffmat2(N,L)
%Derivative matrix of first order
%  
D2 = zeros(N+1,N+1);
for i = 1:N+1
       for j = i+2:2:N+1
           D2(i,j) =(j-1)*((j-1)^2-(i-1)^2);
       end
    end
    D2(1,:) = D2(1,:)/2;
    D2 = D2/(L/2)^2;