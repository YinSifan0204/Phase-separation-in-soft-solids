%%%%%%%%%%%%%%%%%
%   Data analysis module
%           Sifan Yin 
%       08/24/2020
  
% clc;clear all;close all;
A= 10;
r0 = 0.5;
D= 0.01;

gamma =0.2;

load('BlueGrayColormap.mat');
paraname = ['Parameters2D_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat']
load(paraname)
dt = para2D(5);tmax = para2D(6);N =  para2D(7);ist = para2D(8);L = para2D(9);
dataname_all =  ['Data0428_ab01_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat'];
load(dataname_all);
M = length(tsol);index = round(M)
filename1 = ['CH0428_ab01_N',num2str(N),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'_T',num2str(index),'.tif'];
filename2 = ['CH0428_ab01_N',num2str(N),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.gif'];
filename3 = ['CH0428_ab01_N',num2str(N),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'_Fspace',num2str(index),'.tif'];



rmin = min(min(rhosol(:,:)))
rmax = max(max(rhosol(:,:)))
dx = L/N;
x = L/N*[-N/2:N/2]; y = x;  % New meshgrid with Right-most side
[X,Y] = meshgrid(x,y);
          

%% Patch figure faces 
mode = 16;  % How many mesh presented 
inter = N/mode;
amp = 1;  % Amplitude of displacement
face =zeros(mode^2,4);
%face(:,1) = 1:inter: (N+1)^2-(N+1);
%face = [];%zeros(N/inter+1,4)
first_row = 1:inter:N+1-inter;
face_row = first_row';
for j = 1:N/inter-1
  every_row =  first_row'+j*inter*(N+1);
  face_row =[face_row;every_row] ;
end

face(:,1) = face_row;
face(:,2) = face(:,1)+(N+1)*inter;
face(:,3) = face(:,2)+inter;
face(:,4) = face(:,1) +inter;



%% Static images 

figure('color',[1 1 1]);clf;
%index=101;
 rho_temp = reshape(rhosol(index,:),N,N);
 u1_temp = reshape(u1sol(index,:),N,N);
 u2_temp = reshape(u2sol(index,:),N,N);

 rho_temp(:,N+1) = rho_temp(:,1);
 rho_temp(N+1,:) = rho_temp(1,:);
 u1_temp(:,N+1) = u1_temp(:,1);
 u1_temp(N+1,:) = u1_temp(1,:);
 u2_temp(:,N+1) = u2_temp(:,1);
 u2_temp(N+1,:) = u2_temp(1,:);

 %colormap cool
 surf(X,Y,rho_temp,'edgecolor','none')
 %colormap(CustomColormap);
 colormap(BlueGrayColormap)
 axis equal
 hold on
 xlim([-L/2 L/2])
 ylim([-L/2 L/2])
 caxis([rmin rmax])
 title(['Time=',num2str(tsol(index))])
 colorbar
 xlabel x,ylabel y,zlabel \rho,view(0,90);
set(gca,'Fontsize',22,'Fontname','Times New Roman')

 
    Xdef = X+u1_temp*amp;
    Ydef = Y+u2_temp*amp;
    xdef = reshape(Xdef,(N+1)^2,1);
    ydef = reshape(Ydef,(N+1)^2,1);
    zdef = ones((N+1)^2,1);
    vertex = [xdef ydef zdef] ;
   % p  = patch('Faces',face,'Vertices',vertex);
   % set(p,'EdgeColor','red','EdgeAlpha',1,'FaceColor','none','FaceAlpha',1,'LineWidth',0.8);
    axis equal
   print(gcf,'-dtiff',filename1)
% 
% 
%  figure('color',[1 1 1]);clf; 
%  for n = 1: M   %(M-1)/10
%      %if mod(n,1/dt/10)==1
%        %ii = 10*(n-1)+1;
%        ii = n
%        rho_temp = reshape(rhosol(ii,:),N,N);
%        rho_temp(:,N+1) = rho_temp(:,1);
%        rho_temp(N+1,:) = rho_temp(1,:);
%        
%         u1_temp = reshape(u1sol(ii,:),N,N);
%         u2_temp = reshape(u2sol(ii,:),N,N);
%         u1_temp(:,N+1) = u1_temp(:,1);
%         u1_temp(N+1,:) = u1_temp(1,:);
%         u2_temp(:,N+1) = u2_temp(:,1);
%         u2_temp(N+1,:) = u2_temp(1,:);
%        
%        surf(x,y,rho_temp,'edgecolor','none');
%        %colormap(CustomColormap);
%        colormap(BlueGrayColormap)
%        axis equal
%        hold on
% %        rmin = min(rhosol(ii,:));
% %        rmax = max(rhosol(ii,:));
%        xlim([-L/2 L/2])
%        ylim([-L/2 L/2])
%        caxis([rmin rmax])
% %       rhomin=min(min(rhosol)); rhomax=max(max(rhosol));
% %       axis([-pi pi -pi pi rhomin rhomax]);
%        title(['Time=',num2str(tsol(ii))])
%        colorbar
%        xlabel x,ylabel y,zlabel \rho,view(0,90);
%        %end
%        set(gca,'Fontsize',22,'Fontname','Times New Roman')
%         
%         Xdef = X+u1_temp*amp;
%         Ydef = Y+u2_temp*amp;
%         xdef = reshape(Xdef,(N+1)^2,1);
%         ydef = reshape(Ydef,(N+1)^2,1);
%         zdef = ones((N+1)^2,1);
%         vertex = [xdef ydef zdef] ;
%        % p  = patch('Faces',face,'Vertices',vertex);
%        % set(p,'EdgeColor','red','EdgeAlpha',1,'FaceColor','none','FaceAlpha',1,'LineWidth',0.8);
%        drawnow;
%         im = frame2im(getframe(gcf));  % Get current figure
%          [AA, map] = rgb2ind(im,256);
%            if (n==1)
%                imwrite(AA,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.0001);
%            else
%                imwrite(AA,map,filename2,'gif','WriteMode','append','DelayTime',0.0001);
%            end
%            clf;
%  end
%  close all
 
 
 figure('color',[1 1 1])
  rho_temp = reshape(rhosol(index,:),N,N);
  rhoh_real = abs(real(fft2(rho_temp)/(N^2)));rhoh_real(1) = 0;
  rhoh_img = imag(fft2(rho_temp)/(N^2));
  kx = (2*pi/L)*[0:N/2-1 -N/2:-1]; ky = kx;
  [kX,kY] = meshgrid(kx,ky);
%   kxintp =  (2*pi/L)*[0:0.1:N/2-1 -N/2:0.1:-1];
%   kyintp = kxintp;
%   [kXintp,kYintp] = meshgrid(kx,ky);
%   rhoh_intp = interp2(kX,kY,rhoh_real,kXintp,kYintp,'spline')
  surf(kX,kY,rhoh_real,'edgecolor','none')
  %colormap(BlueGrayColormap)
  %colormap(CustomColormap);
  %surf(kXintp,kYintp,rhoh_intp,'edgecolor','none')
   xlim([-N*pi/L/8 N*pi/L/8])
   ylim([-N*pi/L/8 N*pi/L/8])
   colorbar
   xlabel kx,ylabel ky,view(0,90);
   [index1, index2] =   find(rhoh_real==max(max(rhoh_real)))
   kXmax = kX(index1(1),index2(1))
   kYmax = kY(index1(1),index2(1))
   kmax1 = sqrt(kXmax^2+kYmax^2)

%    find(rhoh_real==min(min(rhoh_real)))
   title(['kmax=', num2str(kmax1)]);
   set(gca,'Fontsize',22,'Fontname','Times New Roman')
   print(gcf,'-dtiff',filename3)
 
 