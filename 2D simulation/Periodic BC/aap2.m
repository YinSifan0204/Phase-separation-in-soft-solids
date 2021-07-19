function ph = aap2( uh, vh) %anti-aliased product, 2D 
% in: uh,vh from fft2 with n samples
n = length(uh);
m = n*3/2;
uhp = zeros(m,m);
vhp = zeros(m,m);
ph  = zeros(n,n); 

uhp(1:n/2,1:n/2) = uh(1:n/2,1:n/2);
uhp(1:n/2,m-n/2+1:m) = uh(1:n/2,n/2+1:n) ;
uhp(m-n/2+1:m,1:n/2) = uh(n/2+1:n,1:n/2);
uhp(m-n/2+1:m,m-n/2+1:m) = uh(n/2+1:n,n/2+1:n);

vhp(1:n/2,1:n/2) = vh(1:n/2,1:n/2);
vhp(1:n/2,m-n/2+1:m) = vh(1:n/2,n/2+1:n) ;
vhp(m-n/2+1:m,1:n/2) = vh(n/2+1:n,1:n/2);
vhp(m-n/2+1:m,m-n/2+1:m) = vh(n/2+1:n,n/2+1:n);

up=ifft2(uhp); vp=ifft2(vhp); w=up.*vp; wh=fft2(w);

wh(n/2+1:m-n/2,:) = [];
wh(:,n/2+1:m-n/2) = [];
ph = 1.5*wh;
end