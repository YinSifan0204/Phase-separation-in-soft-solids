   %%%%%%%%%%%%%%%%%
    %
    %        FCT transformation test
    %         x-Fourier    y-Chebyshev
    %            0--------->Y
    %            |
    %            |
    %            |
    %            X
    %
    %            Sifan Yin 
    %         09/04/2020
    %%%%%%%%%%%%%%%%%%
    clc;clear all;close all
    
%% Parameters   

tic;
A= 5;
D= 0.01;
r0 = 0.5;
gamma = 0.1;

Lx = 5; %2*pi;
Ly = 10;
Nx = 4;
Ny = 8;

dataname_all = ['CH2D_0703_FCT_ab1v2001_Lx',num2str(Lx),'_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat'];


dt = 0.005;
tmax =10;
ist = 0.1/dt;



rhomin = 0;
rhomax = 1;


%% Mechanical Parameters
alpha = -1; beta = 1; 
Ga = 0.1;Gp = 1;lambdaa = 0.1;
mu1a = 0;mu1p = 0;
mu2a = 0.1;
m=2; K=2;    
penalty = 0.001;
Gamma = dt;
amp =1;
 
    %% Preparation
    % mesh grids
    dx = Lx/Nx;
    x   =  dx*[-Nx/2:Nx/2-1];
    dtheta = pi/Ny;
    theta = pi:-dtheta:0;
    y = cos(theta)*Ly/2;
    [X,Y] = ndgrid(x,y);   %derivate separately
    
    kx = (2*pi/Lx)*[0:Nx/2-1 -Nx/2:-1];
    kX = kx.'.*ones(1,Ny+1);
    Df1 = diag(1i*kx);
    Df2 = diag(-kx.^2);
    Df4 = diag(kx.^4);

    Dch1 = chebDiffmat(Ny,Ly);
    Dch2 = chebDiffmat2(Ny,Ly);
    Dch3 = Dch1^3;
    Dch4 = Dch2^2;
  
    Ix= eye(Nx,Nx);
    Iy = eye(Ny+1,Ny+1);

    Dx = kron(Iy,Df1);    % Make sure! ref: all_function_test.m
    Dx2 = kron(Iy,Df2);
    Dx4 = kron(Iy,Df4);
    Dy  = kron(Dch1,Ix); 
    Dy2  = kron(Dch2,Ix);
    Dy4  = kron(Dch4,Ix);
    Dxy = Dx*Dy;
    Dxy2 =Dx2*Dy2;
    Dxy4 = Dx4+Dy4+2*Dxy2;
    %Dxy4 = (Dx2+Dy2)^2;
    Iden = eye(Nx*(Ny+1),Nx*(Ny+1));
    PreD_time = toc
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initial conditions--P represents physical space
      tic;
    rho0_P = zeros(Nx,Ny+1);  
    u10_P = zeros(Nx,Ny+1);
    v10_P = zeros(Nx,Ny+1);
    u20_P = zeros(Nx,Ny+1);
    v20_P = zeros(Nx,Ny+1);
    
   nx = 1;   ny = 1;  ampt = 0.01;
  for j1 = 0 : Nx/nx-1 
        for j2 = 0 : Ny/ny-1
            rho0_P(j1*nx+1:(j1+1)*nx, j2*ny+1:(j2+1)*ny)=ampt*(rand(1,1)-0.5)*ones(nx,ny)+r0; 
        end
  end
  
 rho0_P(:,1) = r0*ones(Nx,1);
 rho0_P(:,Ny+1) = rho0_P(:,1);
 
 
% Boundary condition

     V1bc = 0;
     V2bc = 0.01;

    rho_y0 = rho0_P(:,1);
    rho_y0h = fft(rho_y0)/Nx/2;
    rho_yNh = rho_y0h;

   v10_P(:,Ny+1) = V1bc*ones(Nx,1);
   v1_y0 = v10_P(:,1);
   v1_yN = v10_P(:,Ny+1);
   v1_y0h = fft(v1_y0)/Nx/2;
   v1_yNh = fft(v1_yN)/Nx/2;
   
   v20_P(:,Ny+1) = V2bc*ones(Nx,1);
   v2_y0 = v20_P(:,1);
   v2_yN = v20_P(:,Ny+1);
   v2_y0h = fft(v2_y0)/Nx/2;
   v2_yNh = fft(v2_yN)/Nx/2;
   
   IC_time = toc
  

%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions
tic;
n = 0:Ny;
o = ones(1,Ny+1);
left0 = (-o).^n ;      % left0*uch=f(x0)
right0 = o ;             %
left1 = left0*Dch1;    %nleft*uch = f'(x0)
right1 = right0*Dch1;
left2 = left0*Dch2;    %nleft*uch = f'(x0)
right2 = right0*Dch2;
left3 = left0*Dch3 ;   %nleft*uch = f'(x0)
right3 = right0*Dch3;


    BC0 = left0;
    BC1 = right0;
    BC2 = left1;
    BC3 = right1;
    
    Dbc0 = kron(BC0,Ix);
    Dbc1 = kron(BC1,Ix);
    Dbc2 = kron(BC2,Ix);
    Dbc3 = kron(BC3,Ix);
    
   % Dbc_rho = [Dbc0;Dbc1;Dbc2;Dbc3];
    
    Dbc_v1     = [Dbc0;Dbc1];
    Dbc_v2     = [Dbc0;Dbc1];
    ZERO = zeros(2*Nx,Nx*(Ny+1));
    %Dbc_v = [Dbc_v1 ZERO; ZERO Dbc_v2];
    
    rhs_rho = zeros(4*Nx,1);
    rhs_rho(1:2*Nx,1) =[rho_y0h;rho_yNh];
    
    rhs_v1 =zeros(2*Nx,1);
    rhs_v2 = [v2_y0h;v2_yNh];

    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coefficients matrix and verification for vector form --success!
 Lrho = Iden + D*dt*gamma*Dxy4;

% Method 2: impose the BCs replacing the highest order terms
   Lrho(Nx*(Ny-3)+1:Nx*(Ny-2),:) = Dbc0;
   Lrho(Nx*(Ny-2)+1:Nx*(Ny-1),:) = Dbc1;
   Lrho(Nx*(Ny-1)+1:Nx*Ny      ,:) = Dbc2;
   Lrho(Nx*Ny+1     :Nx*(Ny+1),:) = Dbc3;
   tic; Lrho_inv = Lrho\Iden;toc  
  
%   Dbc_rho = [Dbc0;Dbc1;Dbc2;Dbc3];
%   Lrho_new = [zeros(4*Nx,Nx*(Ny+1));Lrho(1:Nx*(Ny-3),:)];spy(Lrho_new)
%   Lrho_new(1:4*Nx,:)=Dbc_rho;
%   tic; Lrhon_inv = Lrho_new\Iden;toc  
%      f_rep(Nx*(Ny-3)+1:Nx*(Ny-2)) = uy0h;
%      f_rep(Nx*(Ny-2)+1:Nx*(Ny-1)) = uyNh;
%       f_rep(Nx*(Ny-1)+1: Nx*Ny        ) = 0;
%      f_rep( Nx*Ny+1      : Nx*(Ny+1)) = 0;



%Lrho_inv = Lrho\eye(Nx*(Ny+1));

% Lrho = [Lrho;Dbc0;Dbc1;Dbc2;Dbc3];
% Lrho_inv = Lrho\eye(Nx*(Ny+5));
% 
 A11 = (1+ dt/Gamma)*Iden  -dt/Gamma*(   (  (2*Gp+1)*dt+( 2*mu1p+1 ) )*Dx2 + (Gp*dt+mu1p)*Dy2  );
 A11(Nx*(Ny-1)+1:Nx*(Ny+1),:) = Dbc_v1;
  
 A12 = -dt/Gamma*( 1 + dt + mu1p + Gp*dt )*Dxy;
 A12(Nx*(Ny-1)+1:Nx*(Ny+1),:)  = ZERO;
 A21 = A12;
 
 A22 = (1+ dt/Gamma)*Iden  -dt/Gamma*(   (  (2*Gp+1)*dt+( 2*mu1p+1 ) )*Dy2 + (Gp*dt+mu1p)*Dx2  );
 A22(Nx*(Ny-1)+1:Nx*(Ny+1),:) = Dbc_v2;
 
 Acmb = [A11 A12 ; A21 A22];
 tic;Acmb_inv = Acmb\eye(2*Nx*(Ny+1));toc
 
% Acmb = [A11 A12 ; A21 A22 ; Dbc_v];
% Acmb_inv = Acmb\eye(2*Nx*(Ny+3));
% tic;[LAcmb,UAcmb] = lu(Acmb);
% tic;Acmb_inv = UAcmb\(  LAcmb \ eye(2*Nx*(Ny+1) ) ); toc
% 


BC_time = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    

    time = 0;
    rho0_FC = FCT(rho0_P);
    u10_FC = FCT(u10_P);  u20_FC = FCT(u20_P);
    v10_FC = FCT(v10_P);  v20_FC = FCT(v10_P);
    
    rho_FC = rho0_FC;
    u1_FC = u10_FC;
    u2_FC = u20_FC;
    v1_FC = v10_FC;
    v2_FC = v20_FC;
    tic
    i = 1; tsol(1)=0;
    while time<tmax
        i = i+1;
        
        for j = 1:ist
            time = time+dt;
            
            %Prep: Matrix form 
            % If Dx,Dy are used, the matrix must be converted to a vector;
            %rho_FCv  = rho_FC(:);
            %u1_FCv = u1_FC(:);
            %u2_FCv = u2_FC(:);
            %v1_FCv = v1_FC(:);
            %v2_FCv = v2_FC(:);
           
            rho_P = iFCT(rho_FC);
           
            
            u1x1_FC = 1i*kX.*u1_FC;         u1x2_FC = u1_FC*Dch1'; 
            u2x1_FC =1i*kX.*u2_FC;          u2x2_FC = u2_FC*Dch1';
                     
            v1x1_FC = 1i*kX.*v1_FC;        v1x2_FC = v1_FC*Dch1';
            v2x1_FC =1i*kX.*v2_FC;         v2x2_FC = v2_FC*Dch1';
           
%            u1x1_FC = reshape(u1x1_FCv,Nx,Ny+1);
%            u2x2_FC = reshape(u2x2_FCv,Nx,Ny+1);
%            u1x2_FC = reshape(u1x2_FCv,Nx,Ny+1);
%            u2x1_FC = reshape(u2x1_FCv,Nx,Ny+1);
           
%             v1x1_FCv = Dx*v1_FCv;
%             v2x2_FCv = Dy*v2_FCv;
%             v1x2_FCv = Dy*v1_FCv;
%             v2x1_FCv = Dx*v2_FCv;
%             
%            v1x1_FC = reshape(v1x1_FCv,Nx,Ny+1);
%            v2x2_FC = reshape(v2x2_FCv,Nx,Ny+1);
%            v1x2_FC = reshape(v1x2_FCv,Nx,Ny+1);
%            v2x1_FC = reshape(v2x1_FCv,Nx,Ny+1);
           
           divU_FC = u1x1_FC + u2x2_FC;
           divV_FC = v1x1_FC + v2x2_FC;
           shearU_FC = u1x2_FC + u2x1_FC;
           shearV_FC = v1x2_FC  +  v2x1_FC;
           
           ds11_FC = 2*(Ga-Gp)*u1x1_FC + 2*(mu1a-mu1p)*v1x1_FC+ (lambdaa-1)*divU_FC+(mu2a-1)*divV_FC;
           se11_FC = aapx(ds11_FC,u1x1_FC);
           
           ds22_FC = 2*(Ga-Gp)*u2x2_FC + 2*(mu1a-mu1p)*v2x2_FC+ (lambdaa-1)*divU_FC+(mu2a-1)*divV_FC;
           se22_FC = aapx(ds22_FC,u2x2_FC);
          
           ds12_FC =    (Ga-Gp)*shearU_FC + (mu1a-mu1p)*shearV_FC;
           se12_FC = aapx(ds12_FC,shearU_FC);
           
           P1 = aapx(ds11_FC,rho_FC);   %P1v = P1(:);
           P2 = aapx(ds12_FC,rho_FC);   %P2v = P2(:);
           Q1 = P2;                                  %Q1v = P2v;
           Q2 = aapx(ds22_FC,rho_FC); %Q2v = Q2(:);
           
         % Equation for rho
         rho_P = iFCT(rho_FC);
        
         
    
         
     

    
            frhop_deno = 1./((1+K*rho_P.^m)).^2; 
            frhop_deno_FC = FCT(frhop_deno);
            frhop_FC = m*aapx(rho_FC,frhop_deno_FC);
            frho_deno = 1./(1+K*rho_P.^m);
            frho_deno_FC = FCT(frho_deno);
            frho_FC = aapx( aapx(rho_FC,rho_FC), frho_deno_FC);
       
            
          % Equations1--rho
          rhov1_FC = aapx(rho_FC,v1_FC);
          %rhov1_FCv = rhov1_FC(:);
          rhov2_FC = aapx(rho_FC,v2_FC);
          %rhov2_FCv = rhov2_FC(:);
          Frho_adv = -dt*( 1i*kX.*rhov1_FC+rhov2_FC*Dch1' );
          
          rho_cub_FC = aapx(aapx(rho_FC,rho_FC),rho_FC);
          Frho_CH= alpha*rho_FC+beta*rho_cub_FC;
          Frho_couple = A*aapx(frhop_FC,divU_FC);
          Frho_mech = 0.5*(se11_FC+se12_FC+se22_FC);
          Frho_diff_FC = Frho_CH + Frho_couple+Frho_mech;
          %Frho_diff_v = Frho_diff_FC(:);
          Frho_diff1 = -dt*D*kX.^2.*Frho_diff_FC;
          Frho_diff2 =  dt*D*Frho_diff_FC*Dch2';
          Frho_diffm= Frho_diff1 + Frho_diff2;
          %Frho_diffv = Frho_diffm(:);
          Frho_m = rho_FC + Frho_adv + Frho_diffm;
          Frho_v = Frho_m(:);
          Frho_v(Nx*(Ny-3)+1:Nx*(Ny+1)) = rhs_rho;
          %Frho_diff_v = dt*D*(Dx2+Dy2)*Frho_diff_v;
          %Frho_v = rho_FCv + Frho_adv_v + Frho_diffv;
          %Frho_aug = [Frho_m(:);rhs_rho];
                 
         % rho_FCv = Lrho_inv*Frho_aug;
          rho_FCv = Lrho_inv*Frho_v;
%           Dbc0*rho_FCv
%           Dbc1*rho_FCv
%           Dbc2*rho_FCv
%           Dbc3*rho_FCv
%           
          
          rho_FC = reshape(rho_FCv,Nx,Ny+1); 
         % rho_FC(1,1) = r0/2;    % Attention!!!
          rho_P = iFCT(rho_FC);
%           rho_P(:,1)
%           rho_P(:,Ny+1)
%           rho0_P(:,1) = r0*ones(Nx,1);
%           rho0_P(:,Ny+1) = rho0_P(:,1);
          
          index_rhomin = find(rho_P<=rhomin);
          rho_P(index_rhomin) = rhomin+penalty;
          index_rhomax = find(rho_P>=rhomax);
          rho_P(index_rhomax) = rhomax - penalty;
          rho_FC = FCT(rho_P);
          rho_FC(1,1) = r0/2;
          rho_P = iFCT(rho_FC);
         
           % Equations2--v1 v2
           PA = A*frho_FC;
%        PAv = PA(:);
%        b1u = dt/Gamma*(  ((2*Gp+1)*Dx2+Gp*Dy2)*u1_FCv+(Gp+1)*Dxy*u2_FCv  );
%        b2u = dt/Gamma*(  ((2*Gp+1)*Dy2+Gp*Dx2)*u2_FCv+(Gp+1)*Dxy*u1_FCv );
            b1u = dt/Gamma*(  -(2*Gp+1)*kX.^2.*u1_FC+ Gp*u1_FC*Dch2'    +(Gp+1)*1i*kX.*u2_FC*Dch1'  );
            b2u = dt/Gamma*(   (2*Gp+1)*u2_FC*Dch2'  -  Gp*kX.^2.*u2_FC   +(Gp+1)*1i*kX.*u1_FC*Dch1'  );
            b1 = v1_FC + dt/Gamma*( 1i*kX.*P1  + P2*Dch1' + 1i*kX.*PA  )+  b1u;
            b2 = v2_FC + dt/Gamma*( 1i*kX.*Q1 + Q2*Dch1'+  PA*Dch1'  ) + b2u;
%            b1 = v1_FCv+dt/Gamma*(Dx*P1v+Dy*P2v+Dx*PAv)+b1u;
%            b2 = v2_FCv+dt/Gamma*(Dx*Q1v+Dy*Q2v+Dy*PAv)+b2u;
           B1 = b1(:); B1(Nx*(Ny-1)+1:Nx*(Ny+1)) = rhs_v1;
           B2 = b2(:); B2(Nx*(Ny-1)+1:Nx*(Ny+1)) = rhs_v2;
           Bcmb = [B1;B2];
           
           vsol = Acmb_inv*Bcmb;     %toc    %% could be improved!
           v1_FCv = vsol(1:Nx*(Ny+1));
           v2_FCv = vsol(Nx*(Ny+1)+1:end);
           v1_FC = reshape(v1_FCv,Nx,Ny+1); %v1_FC(1,1) = 0;
           v2_FC = reshape(v2_FCv,Nx,Ny+1); %v2_FC(1,1) = 0;
           
            u1_FC = u1_FC + dt*v1_FC; 
            u2_FC = u2_FC + dt*v2_FC; %u1h(1,1) = 0;
            u1_P = iFCT(u1_FC);
            u2_P = iFCT(u2_FC);   %u2h(1,1) = 0;
            v1_P = iFCT(v1_FC);
            v2_P = iFCT(v2_FC);
            v_mag = sqrt(v1_P.^2+v2_P.^2);
            vmax = max(max(v_mag));
            cfl = vmax*dt/(1.41*dx);
          if cfl>0.5
               cfl
               time
               dt = dt/2;
             break
          end
        end
       time
       tsol(i) = time;
       rhosol(i,:) = reshape(rho_P,1,Nx*(Ny+1));
       u1sol(i,:) = reshape(u1_P,1,Nx*(Ny+1));
       u2sol(i,:) = reshape(u2_P,1,Nx*(Ny+1));
       %Iteration_time = toc
    end
    toc;
rmin = min(min(rhosol(:,:)));
rmax = max(max(rhosol(:,:)))    ;


tsol = tsol';
save(dataname_all,'tsol','rhosol','u1sol','u2sol')
toc




%%%%%%%%%%%%%%%%%%%%%
%   Post-processing
% load(dataname_all)
% M = size(rhosol,1);
% 
% 
%  %% Patch figure faces 
%  
% mode_x = 16;  % How many mesh presented 
% mode_y = 32;
% 
% inter_x = Nx/mode_x;
% inter_y = Ny/mode_y;
% amp = 1;  % Amplitude of displacement
% 
% face =zeros(mode_x*mode_y,4);
% first_row = 1:inter_y:Ny+1-inter_y;
% face_row = first_row';
% for j = 1:Nx/inter_x-1
%   every_row =  first_row'+j*inter_x*(Ny+1);
%   face_row =[face_row;every_row] ;
% end
% 
% face(:,1) = face_row;
% face(:,2) = face(:,1) + (Ny+1)*inter_x;
% face(:,3) = face(:,2) + inter_y;
% face(:,4) = face(:,1) + inter_y;
% 
% 
% 
%  index = round(M)-1;
% if (BC_int==1)
%    filename1 = ['CH2D_FCTperi_ab1','_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'_T',num2str(index),'.tif'];
%    filename2 = ['CH2D_FCTperi_ab1','_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.gif'];
% else
%    filename1 = ['CH2D_FCTdn_ab1v2','_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'_T',num2str(index),'.tif'];
%    filename2 = ['CH2D_FCTdn_ab1v2','_A',num2str(A),'_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.gif'];
% end  
% %% Static images
%  figure(1);clf;
%   
%  rho_temp = reshape(rhosol(index,:),Nx,Ny+1);
%  u1_temp = reshape(u1sol(index,:),Nx,Ny+1);
%  u2_temp = reshape(u2sol(index,:),Nx,Ny+1);
%  
%  rho_temp(Nx+1,:) = rho_temp(1,:);
%  u1_temp(Nx+1,:) = u1_temp(1,:);
%  u2_temp(Nx+1,:) = u2_temp(1,:);
%  
% rmin = min(min(rho_temp));
% rmax = max(max(rho_temp));
%  %colormap cool
%  x = Lx/Nx*[-Nx/2:Nx/2]; % New meshgrid with Right-most side
%  [X,Y] =  meshgrid(x,y);
%  xx = -Lx/2:0.01:Lx/2; yy = -Ly/2:0.01:Ly/2;
%  [XX,YY] = meshgrid(xx,yy);
% 
% rho_intp = interp2(X,Y,rho_temp',XX,YY,'spline');
% surf(XX,YY,rho_intp,'edgecolor','none')
%  
% %   u1_intp = interp2(X,Y,u1_temp',XX,YY,'spline');
% %   surf(XX,YY,u1_intp,'edgecolor','none')
% %  u2_intp = interp2(X,Y,u2_temp',XX,YY,'spline');
% %  surf(XX,YY,u2_intp,'edgecolor','none')
% %[X,Y] =  ndgrid(x,y);
% %surf(X,Y,rho_temp,'edgecolor','none')
%  axis equal
%  hold on
%  xlim([-Lx/2 Lx/2])
%  ylim([-Ly/2 Ly/2])
%  caxis([rmin rmax])
%  title(['Time=',num2str(tsol(index))])
%  colorbar
%  xlabel x,ylabel y, zlabel \rho,
%  view(90,90);
% set(gca,'Fontsize',22,'Fontname','Times New Roman')
% 
%     Xdef = X+u1_temp'*amp;
%     Ydef = Y+u2_temp'*amp;
%      xp = Lx/Nx*[-Nx/2:Nx/2];
%      yp =Ly/Ny*[-Ny/2:Ny/2];
%      [Xp,Yp] = meshgrid(xp,yp);
%      Xdef = interp2(X,Y,Xdef,Xp,Yp,'spline');
%      Ydef = interp2(X,Y,Ydef,Xp,Yp,'spline');
%     xdef = reshape(Xdef,(Nx+1)*(Ny+1),1);
%     ydef = reshape(Ydef,(Nx+1)*(Ny+1),1);
%     zdef = ones((Nx+1)*(Ny+1),1);
%     vertex = [xdef ydef zdef] ;
%     p  = patch('Faces',face,'Vertices',vertex);
%    set(p,'EdgeColor','red','EdgeAlpha',1,'FaceColor','none','FaceAlpha',1,'LineWidth',0.5);
%    set(gcf,'position',[500 250 1000 500])
%    print(gcf,'-dtiff',filename1)  
%     
%    
%    figure(2);
%      for n = 1: M   %(M-1)/10
%        if mod(n,10)==1
%        ii = n
%        %ii = n
%        rho_temp = reshape(rhosol(ii,:),Nx,Ny+1);
%        u1_temp = reshape(u1sol(ii,:),Nx,Ny+1);
%        u2_temp = reshape(u2sol(ii,:),Nx,Ny+1);
%  
%       rho_temp(Nx+1,:) = rho_temp(1,:);
%       u1_temp(Nx+1,:) = u1_temp(1,:);
%       u2_temp(Nx+1,:) = u2_temp(1,:);
%        
%       rho_intp = interp2(X,Y,rho_temp',XX,YY,'spline');
%       surf(XX,YY,rho_intp,'edgecolor','none')
%       axis equal
%       hold on
%       xlim([-Lx/2-1 Lx/2+1])
%       ylim([-Ly/2-1 Ly/2+1])
%       caxis([rmin rmax])
%      title(['Time=',num2str(tsol(ii))])
%      colorbar
%     xlabel x,ylabel y, zlabel \rho,
%     view(90,90);
%      set(gca,'Fontsize',22,'Fontname','Times New Roman')
% 
%     Xdef = X+u1_temp'*amp;
%     Ydef = Y+u2_temp'*amp;
%      xp = Lx/Nx*[-Nx/2:Nx/2];
%      yp =Ly/Ny*[-Ny/2:Ny/2];
%      [Xp,Yp] = meshgrid(xp,yp);
%      Xdef = interp2(X,Y,Xdef,Xp,Yp,'spline');
%      Ydef = interp2(X,Y,Ydef,Xp,Yp,'spline');
%     xdef = reshape(Xdef,(Nx+1)*(Ny+1),1);
%     ydef = reshape(Ydef,(Nx+1)*(Ny+1),1);
%     zdef = ones((Nx+1)*(Ny+1),1);
%     vertex = [xdef ydef zdef] ;
%     p  = patch('Faces',face,'Vertices',vertex);
%     set(p,'EdgeColor','red','EdgeAlpha',1,'FaceColor','none','FaceAlpha',1,'LineWidth',0.5);
%     set(gcf,'position',[500 250 1000 500])
%        drawnow;
%         im = frame2im(getframe(gcf));  % Get current figure
%          [AA, map] = rgb2ind(im,256);
%            if (n==1)
%                imwrite(AA,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.0001);
%            else
%                imwrite(AA,map,filename2,'gif','WriteMode','append','DelayTime',0.0001);
%            end
%            clf;
%        end
%      end
%  close all
    
    
    
    
    
    
    
    
    
    
    