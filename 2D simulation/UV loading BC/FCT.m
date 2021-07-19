function u_fch = FCT(u)
% x--Fourier  y--Chebyshev
%
Nx = size(u,1);
Ny = size(u,2)-1;
u_f = fft(u,[],1)/Nx;  % FCT and iFCT use  ' /Nx  ,*Nx'
u_fch = chyfft(u_f);

end

