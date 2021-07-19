function [nleft,nright] = NeumannBC(N,L)
% Neumann Boundary
 n = 0:N;
 o = ones(1,N+1);
 dleft = (-o).^n;  % dleft*uch=f(x0)
 dright = o ;  %dright*uch = f(xN)
 D1 = chebDiffmat(N,L);
 nleft = dleft*D1;    %nleft*uch = f'(x0)
 nright = dright*D1;
end

