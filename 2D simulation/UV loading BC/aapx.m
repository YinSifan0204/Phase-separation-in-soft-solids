function  wh= aapx(uh,vh)
% for FCT method calculate anti-aliasing product 
%  padding every column
Nx = size(uh,1);
Ny = size(uh,2)-1;
M = Nx*3/2;
wh = zeros(Nx,Ny+1);

uh_pad = zeros(M,Ny+1);
uh_pad(1:Nx/2,:) = uh(1:Nx/2,:);
uh_pad(M-Nx/2+1:M,:) = uh(Nx/2+1:Nx,:);
up_pad = iFCT(uh_pad);

vh_pad = zeros(M,Ny+1);
vh_pad(1:Nx/2,:) = vh(1:Nx/2,:);
vh_pad(M-Nx/2+1:M,:) = vh(Nx/2+1:Nx,:);
vp_pad = iFCT(vh_pad);

wp_pad = up_pad.*vp_pad;
wh_pad = FCT(wp_pad);
wh(1:Nx/2,:) = wh_pad(1:Nx/2,:);
wh(Nx/2+1:Nx,:) = wh_pad(M-Nx/2+1:M,:);






end

