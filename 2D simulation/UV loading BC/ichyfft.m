function u_p = ichyfft(uch)
% From Chebyshev space to physical space
Nx = size(uch,1);
Ny = size(uch,2)-1;

M = zeros(Nx,2*Ny);

uch(:,1) = uch(:,1)*2;
uch(:,end) = uch(:,end)*2;

M(:,1:Ny+1) = uch;
M(:,Ny+2:end) = uch(:, end-1:-1:2);
v = ifft(M,[],2);
u_p = v(:,Ny+1:-1:1)*(2*Ny);

end

