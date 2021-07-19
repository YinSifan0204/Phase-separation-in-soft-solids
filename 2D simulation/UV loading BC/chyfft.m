function uch= chyfft(u)
% Chebyshev expansion truncation coefficients
% For 2D matrix --y chebyshev interpolation
Nx = size(u,1);
Ny = size(u,2)-1;

v = zeros(Nx,2*Ny);
v(:,1:Ny+1) =u(:,end:-1:1);
v(:,Ny+2:end) = u(:,2:Ny);

a = fft(v,[],2)/(2*Ny);   % fft  /Ny  ifft *Ny

a = a(:,1:Ny+1);
a(:,1) = a(:,1)/2;
a(:,Ny+1) = a(:,Ny+1)/2;

uch = a;

% if (mod(Nx,2)==0)
%     uch(Nx/2+1,:)=0;
% end

end

