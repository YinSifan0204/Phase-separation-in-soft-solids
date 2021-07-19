%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              CH-mechanical coupled equation
%     Spectral method with semi-implicit time marching
%                Dealiasing--padding aap2.m
%                  debug version
%                       2D model
%                       Sifan Yin
%                     
%                   2021/4/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;

%% para = [A, D, r0, gamma,dt,tmax,N];
A= 10;
r0 = 0.5;

D= 0.01;
gamma =0.01;

paraname = ['Parameters2D_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat']
load(paraname)
alpha = -0.1; beta = 0.1; 
dt = para2D(5);tmax = para2D(6);N =  para2D(7);ist = para2D(8);L = para2D(9);
dataname_all =  ['Data0428_ab01_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat'];


Lx = L;
Ly = L;
dx = Lx/N;
x = Lx/N*[-N/2:N/2-1]; 
y = Ly/N*[-N/2:N/2-1]; 
kx = (2*pi/Lx)*[0:N/2-1 -N/2:-1]; 
ky = (2*pi/Ly)*[0:N/2-1 -N/2:-1]; 
[X,Y] = meshgrid(x,y);
[kX,kY] = meshgrid(kx,ky);
k2 = kX.^2+kY.^2;
k4 = k2.^2;

%% Parameters

m=2; K=2;
Ga = 0.1;Gp = 1;lambdaa = 0.1;
%mu1a = 0.1;mu1p = 1;mu2a = 0.1;
%20200816
mu1a = 0;mu1p = 0;
mu2a = 0.1;
rhomin = 0;
rhomax = 1;        % D=0.2时改变为 0.3~0.7
penalty = 0.001;

% Time and time step
%dt = 0.01; 
%Gamma = dt;
% 2020.8.11 modify Gamma a little bit bigger
Gamma = dt;



%% Algorithm:
%  nonlinear terms
%   --CH  aap
%  --others real(ifft()) as usual
%  Linear terms -- to the left side 

cofrho = 1./(1+dt*D*gamma*k4);
a11 = 1+ dt/Gamma*(   1+kX.^2*(  (2*Gp+1)*dt+( 2*mu1p+1 ) ) + kY.^2*(Gp*dt+mu1p)  );
a12 = dt/Gamma*( 1 + dt + mu1p + Gp*dt )*kX.*kY;
a21 = a12;
a22 = 1+ dt/Gamma*(   1+kY.^2*(  (2*Gp+1)*dt+( 2*mu1p+1 ) ) + kX.^2*(Gp*dt+mu1p)  );
deno = a11.*a22 - a12.*a21;
%cofv = 1./(1+dt/Gamma+dt*(1+dt)/Gamma*k2);
nb = 8;
rho0   =  zeros(N,N);  ampt = 0.1;
u10    =  zeros(N,N); u20 = zeros(N,N);
v10    =  zeros(N,N); v20  =zeros(N,N);  

% %Random initial distribution
 for j1 = 0 : N/nb-1 
    for j2 = 0 : N/nb-1
        rho0(j1*nb+1:(j1+1)*nb, j2*nb+1:(j2+1)*nb)=ampt*(rand(1,1)-0.5)*ones(nb,nb)+r0; 
        %  rho0(j1*nb+1:(j1+1)*nb, j2*nb+1:(j2+1)*nb) = 0.5*(rand(1,1)-0.5)*ones(nb,nb);
    end
 end
%rho0 = u0;
% concentrated initial distribution
%  rho0 = (r0-ampt)*ones(N,N)   ;
%  rho0(7*N/16+1 : 9*N/16, 7*N/16+1 : 9*N/16)= r0 + ampt;

 t = 0;
 %tv = 0 : dt*ist : tmax;
 %tn = length(tv);
 rhoh = fft2(rho0);
 %rhoh(1,1)=N^2*r0;   %保证平均值不变,二维时为N^2,一维时为N
 u1h = fft2(u10); u2h = fft2(u20);
 v1h = fft2(v10);  v2h  =fft2(v20);
 
%%  General method to convert matrix to vector
%   matrix = rand(4,4)
%   vector = matrix(:)
%   Mrecover = reshape(vector,4,4)
%   matrix-Mrecover
tic
time = 0;
tsol(1) = 0; i=1;
rhosol(1,:)  = reshape(rho0,1,N^2);
u1sol(1,:)   = reshape(u10,1,N^2);
u2sol(1,:)   = reshape(u20,1,N^2);
v1sol(1,:)   = reshape(v10,1,N^2);
v2sol(1,:)   = reshape(v20,1,N^2);    % Stored in row
 while time<tmax
     i = i+1;
     for j = 1:ist
          time = time+dt
         % 1. Equation of rhoh
          Frho_adv = - dt*1i*kX.*aap2(rhoh,v1h) - dt*1i*kY.*aap2(rhoh, v2h);
         
         % Miekte
         %frhop_demo = 1./((r0^m+K*rho.^m)).^2; 
         %Hill's function
         rho = real(ifft2(rhoh));
         frhop_deno = 1./((1+K*rho.^m)).^2; 
         frhoph_deno = fft2(frhop_deno);
         frhoph = m*aap2(rhoh,frhoph_deno);
         frho_deno = 1./(1+K*rho.^m);
         frhoh_deno = fft2(frho_deno);
         frhoh = aap2( aap2(rhoh,rhoh), frhoh_deno);
         
         u1x1h = 1i*kX.*u1h;
         u2x2h = 1i*kY.*u2h;
         u1x2h = 1i*kY.*u1h;
         u2x1h = 1i*kX.*u2h;
         
         v1x1h = 1i*kX.*v1h;
         v2x2h = 1i*kY.*v2h;
         v1x2h = 1i*kY.*v1h;
         v2x1h = 1i*kX.*u2h;
         
         divUh = u1x1h + u2x2h;
         divVh = v1x1h + v2x2h;
         shearUh = u1x2h + u2x1h;
         shearVh = v1x2h + v2x1h;
         sh11 = 2*(Ga-Gp)*u1x1h + 2*(mu1a-mu1p)*v1x1h+(lambdaa-1)*divUh+(mu2a-1)*divVh;
         seh11 = aap2(sh11,u1x1h);
         sh22 = 2*(Ga-Gp)*u2x2h+ 2*(mu1a-mu1p)*v2x2h+(lambdaa-1)*divUh+(mu2a-1)*divVh;
         seh22 = aap2(sh22,u2x2h);
         sh12 =    (Ga-Gp)*shearUh + (mu1a-mu1p)*shearVh;
         seh12 = aap2(sh12,shearUh);
         
%          seh11 = 2*(Ga-Gp)*aap2(u1x1h,u1x1h) + 2*(mu1a-mu1p)*aap2(v1x1h,u1x1h)+(lambdaa-1)*aap2(divUh,u1x1h)+(mu2a-1)*aap2(divVh,u1x1h);
%          seh22 = 2*(Ga-Gp)*aap2(u2x2h,u2x2h) + 2*(mu1a-mu1p)*aap2(v2x2h,u2x2h)+(lambdaa-1)*aap2(divUh,u2x2h)+(mu2a-1)*aap2(divVh,u2x2h);
%          seh12 =    (Ga-Gp)*aap2(shearUh,shearUh) + (mu1a-mu1p)*aap2(shearVh,shearUh);
         
         Frho_CH = alpha*rhoh +beta*aap2(aap2(rhoh,rhoh),rhoh);
         Frho_couple = A*aap2(frhoph,divUh);
         Frho_mech = 0.5*(seh11+seh22+seh12);
         Frho_diff = - dt*D*k2.*(Frho_CH + Frho_couple + Frho_mech);
         
         
     
       
         
        % 2. Equations of v1h, v2h 
%          PQ =  (lambdaa-1)*aap2(rhoh,divUh) + (mu2a-1)*aap2(rhoh,divVh);
%          P1 =  2*(Ga-Gp)*aap2(rhoh,u1x1h) + 2*(mu1a-mu1p)*aap2(rhoh,v1x1h) +PQ;
%          P2 =    (Ga-Gp)*aap2(rhoh,shearUh)+ (mu1a-mu1p)*aap2(rhoh,shearVh);
%          Q1 = P2;
%          Q2 = 2*(Ga-Gp)*aap2(rhoh,u2x2h) + 2*(mu1a-mu1p)*aap2(rhoh,v2x2h) + PQ;
          P1 = aap2(sh11,rhoh);
          P2 = aap2(sh12,rhoh);
          Q1 = P2;
          Q2 = aap2(sh22,rhoh);
        
         PA = A*frhoh;
         b1h_u = dt/Gamma*(  ( (2*Gp+1)*kX.^2 + Gp*kY.^2 ).*u1h + (Gp+1)*kX.*kY.*u2h ); 
         b2h_u = dt/Gamma*(  ( (2*Gp+1)*kY.^2 + Gp*kX.^2  ).*u2h+ (Gp+1)*kX.*kY.*u1h  ); 
         b1h = v1h + dt/Gamma*( 1i*kX.*P1 + 1i*kY.*P2  +  1i*kX.*PA ) - b1h_u;
         b2h = v2h + dt/Gamma*( 1i*kX.*Q1 + 1i*kY.*Q2 +  1i*kY.*PA) - b2h_u; 
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Iterations Equations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               rhoh =  cofrho.*( rhoh  + Frho_diff );
          rhoh =  cofrho.*( rhoh + Frho_adv + Frho_diff );
%          rhoh(1,1) = N^2*r0;   % average value is restricted
          rho_real = real(ifft2(rhoh));
         index_rhomin = find(rho_real<=rhomin);
         rho_real(index_rhomin) = rhomin+penalty;
         index_rhomax = find(rho_real>=rhomax);
         rho_real(index_rhomax) = rhomax - penalty;
         rhoh = fft2(rho_real);
         rhoh(1,1) = N^2*r0;
         rho_real = real(ifft2(rhoh));
     
     
         v1h =  (a22.*b1h - a12.*b2h)./deno;% v1h(1,1) = 0; 
         v2h =  (a11.*b2h - a21.*b1h)./deno;% v2h(1,1) = 0;
         u1h = u1h + dt*v1h; %u1h(1,1) = 0;
         u2h = u2h + dt*v2h; %u2h(1,1) = 0;
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          u1h = fft2(u10); u2h = fft2(u20);
%          v1h = fft2(v10);  v2h  =fft2(v20);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         u1_real = real(ifft2(u1h));
         u2_real = real(ifft2(u2h));
         v1_real = real(ifft2(v1h));
         v2_real = real(ifft2(v2h));
         v_mag = sqrt(v1_real.^2+v2_real.^2);
         vmax = max(max(v_mag));
         cfl = vmax*dt/dx;
         if cfl>0.4
             cfl
             time
             dt = dt/2;
             break
         end
    
     end
     time
     tsol(i)       = time;
     rhosol(i,:) = reshape(rho_real,1,N^2);
     u1sol(i,:)  = reshape(u1_real,1,N^2);
     u2sol(i,:)  = reshape(u2_real,1,N^2);
     v1sol(i,:)  = reshape(v1_real,1,N^2);
     v2sol(i,:)  = reshape(v2_real,1,N^2);
%      figure(1)
%        surf(x,y,rho_real,'edgecolor','none')
%        xlim([-pi pi-dx])
%        ylim([-pi pi-dx])
%        caxis([0 1])
%        title(['Time=',num2str(time)])
%        colorbar
%        xlabel x,ylabel y,zlabel \rho,view(0,90);
%        set(gca,'Fontsize',14,'Fontname','Times New Roman')
      
 end


tsol = tsol';
save(dataname_all,'tsol','rhosol','u1sol','u2sol')
toc
 
% M = length(tsol);
% index = round(2*(M-1)/3);
%  figure(1);clf;
%  filename1 = ['Hill2D_N',num2str(N),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'-T',num2str(index),'.tif'];
%  surf(x,y,reshape(rhosol(index,:),N,N),'edgecolor','none')
%  xlim([-L/2 L/2-dx])
%  ylim([-L/2 L/2-dx])
%  caxis([0 1])
%  title(['Time=',num2str(tsol(index))])
%  colorbar
%  xlabel x,ylabel y,zlabel \rho,view(0,90);
% set(gca,'Fontsize',14,'Fontname','Times New Roman')
% print(gcf,'-dtiff',filename1)
%  
% 
% 
% 
%  figure(2);clf; 
%  filename = ['Hill2Drho_N',num2str(N),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'.gif'];
%  for n = 1:(M-1)/10
%      %if mod(n,1/dt/10)==1
%        ii = 10*(n-1)+1
%        surf(x,y,reshape(rhosol(ii,:),N,N),'edgecolor','none')
%        rmin = min(rhosol(ii,:));
%        rmax = max(rhosol(ii,:))
%        xlim([-pi pi-dx])
%        ylim([-pi pi-dx])
%        caxis([rmin rmax])
% %       rhomin=min(min(rhosol)); rhomax=max(max(rhosol));
% %       axis([-pi pi -pi pi rhomin rhomax]);
%        title(['Time=',num2str(tsol(ii))])
%        colorbar
%        xlabel x,ylabel y,zlabel \rho,view(0,90);
%        %end
%        set(gca,'Fontsize',14,'Fontname','Times New Roman')
%        drawnow;
%         im = frame2im(getframe(gcf));  % Get current figure
%          [AA, map] = rgb2ind(im,256);
%            if (n==1)
%                imwrite(AA,map,filename,'gif','LoopCount',Inf,'DelayTime',0.0001);
%            else
%                imwrite(AA,map,filename,'gif','WriteMode','append','DelayTime',0.0001);
%            end
%  end
%  
%  figure(3);clf;
%  M = length(tsol)
%  filename = ['CHmech2D_u1sol_N',num2str(N),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'.gif'];
%  for n = 1:M
%      %if mod(n,1/dt/10)==1
%        surf(x,y,reshape(u1sol(n,:),N,N),'edgecolor','none') 
%        u1min=min(min(u1sol)); u1max=max(max(u1sol));
%        xlim([-pi pi-dx])
%        ylim([-pi pi-dx])
%        caxis([u1min u1max])
%       
% %       axis([-pi pi -pi pi rhomin rhomax]);
%        title(['Time=',num2str(tsol(n))])
%        colorbar
%        xlabel x,ylabel y,zlabel \rho,view(0,90);
%        %end
%        set(gca,'Fontsize',14,'Fontname','Times New Roman')
%        drawnow;
%         im = frame2im(getframe(gcf));  % Get current figure
%          [AA, map] = rgb2ind(im,256);
%            if (n==1)
%                imwrite(AA,map,filename,'gif','LoopCount',Inf,'DelayTime',0.0001);
%            else
%                imwrite(AA,map,filename,'gif','WriteMode','append','DelayTime',0.0001);
%            end
%  end
 










