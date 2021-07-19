function u_re = iFCT(u_fch)
% x--Fourier  y--Chebyshev

Nx = size(u_fch,1);
Ny = size(u_fch,2)-1;
u_fp = ichyfft(u_fch);
u_re = real(ifft(u_fp,[],1))*Nx;

end

