%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CH-mechanical coupled equation
%     Spectral method with semi-implicit time marching
%     Dealiasing--padding aap.m
%     Hill's function = A*rho^m/(1+K*rho^m)
%                 Sifan Yin
%               2020/02/17-2/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;
% Space 
%% para = [A, D, r0, gamma,dt,tmax,N];

A=4
r0 = 0.6
D= 0.01;
paraname = ['Parameters1D0817_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'.mat']
load(paraname)
gamma = para(4); 
dt = para(5);
tmax = para(6);
N = para(7);


L = 10;
%N =2048;   % small D 
dx = L/N;
x = L/N*[-N/2:N/2-1];
kX = (2*pi)/L*[0:N/2-1 -N/2:-1];
k2 = kX.^2;
k4 = k2.^2;

%% Parameters
%A= 5;
%D= 0.0001;
%r0 = 0.5;
alpha = -1; beta = 1; gamma = 0.1;
Ea = 0.1; zetaa = 0.1;
m=2; K=2;
rhomin = 0;
rhomax = 1;
penalty = 0.001;

%dt = 0.0001;
Gamma = 1;
ist = para(8);  % Output interval
%tmax = 20;



%% Algorithm:
%  nonlinear terms
%   --CH  aap
%  --others real(ifft()) as usual
%  Linear terms -- to the left side 

cofrho = 1./(1+dt*D*gamma*k4);
cofv = 1./(1+dt/Gamma+dt*(1+dt)/Gamma*k2);
nb = 1;
rho0 = zeros(1,N);  ampt = 0.1;
u0  = zeros(1,N);
v0  = zeros(1,N);
% sigma = 4;
% rho0 = exp(-0.5*(x.^2)/sigma^2);
 for j1 = 0 : N/nb-1
          rho0(j1*nb+1:(j1+1)*nb) = (r0+ampt*(rand(1,1)-0.5))*ones(1,nb);
 end

% rho0 = (r0-ampt)*ones(1,N)   ;
% rho0(126*N/256+1 : 130*N/256)=r0+ampt;


 t = 0;
 tv = 0 : dt*ist : tmax;
 tn = length(tv);
 rhosol = zeros(tn, N);
 usol = zeros(tn,N);
 vsol = zeros(tn,N);
 rhosol(1,:) = rho0;
 usol(1,:) = u0;
 vsol(1,:) = v0;
 time = 0;
 rhoh = fft(rho0); rhoh(1)=N*r0;
 uh = fft(u0);
 vh = fft(v0);
 tic
for i = 1:tn
    for j = 1:ist
         % input--real space function
         time = time + dt;
         Frho1 = aap(rhoh,vh);
         Frho2_CH = alpha*rhoh +beta*aap(aap(rhoh,rhoh),rhoh);
         

%  Original Hill's function
         rho = real(ifft(rhoh));
         frhop_deno = 1./((1+K*rho.^m).^2);
         frhoph_deno = fft(frhop_deno);
         frhoph = m*A*aap(rhoh,frhoph_deno);  %m=2
         
         uxh = 1i*kX.*uh;
         vxh = 1i*kX.*vh;
%      ux = real(ifft(uxh));
%      vx = real(ifft(vxh));
         Frho2_couple = aap(frhoph, uxh);
         Frho2_ux2 = 0.5*(Ea-1)*aap(uxh,uxh);
         Frho2_uxvx = 0.5*(zetaa-1)*aap(uxh,vxh);
         Frho2 = Frho2_CH+Frho2_couple+Frho2_ux2+Frho2_uxvx;
         
         %Ref: Weber-2018
         %Fv_frho = A*(  rhoh_err  -a*aap( aap(rhoh_err,rhoh_err), rhoh_err) );
         
         % Ref: Mietke-2019
         frho_deno = 1./(1+K*rho.^m);
         frhoh_deno = fft(frho_deno);
         Fv_frho = A*aap( aap(rhoh,rhoh), frhoh_deno);
         Fv_ux = (Ea-1)*aap(rhoh,uxh);
         Fv_vx = (zetaa-1)*aap(rhoh,vxh);
         Fv = Fv_frho + Fv_ux + Fv_vx;

         rhoh = cofrho.*( rhoh - 1i*dt*kX.*Frho1 - D*dt*k2.*Frho2 );  rhoh(1) = N*r0;
         rho_real = real(ifft(rhoh));
         index_rhomin = find(rho_real<rhomin);
         rho_real(index_rhomin) = rhomin + penalty;
         index_rhomax = find(rho_real>rhomax);
         rho_real(index_rhomax) = rhomax - penalty;
         rhoh = fft(rho_real);rhoh(1) = N*r0;
           
         vh     = cofv.*( vh - dt/Gamma*k2.*uh + dt/Gamma*1i*kX.*Fv); vh(1)=0;
         uh = uh + dt*vh; uh(1)=0;
         u_real = real(ifft(uh));
         v_real = real(ifft(vh));        
         vmax = max(v_real);
         cfl = vmax*dt/dx;  % verify if cfl<0.5
         if cfl>0.4
             cfl
             time
             break
         end
         
     end
     time
      rhosol(i,:) = rho_real;
      usol(i,:) = u_real;
      vsol(i,:) = v_real;
 end
toc

dataname =  ['Data1D0817_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat'];
save(dataname,'x','tv','rhosol','usol','vsol');


%  M = length(tv);
% for i = 1:M
%         mean_rho(i) = mean(rhosol(i,:));
%         mean_moment(i) = mean(rhosol(i,:).*x);
%         mean_square(i)= mean(rhosol(i,:).*(x.^2));
% end
% 
%  
% figure(1);clf;
% filename1 = ['CH_July22_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),...
%                        '_gamma',num2str(gamma),'_T',num2str(tmax),'.tif'];
%  x(N+1) = L/2; 
%  [X,T] = meshgrid(x,tv);
%  rhosol(:,N+1) = rhosol(:,1);
% surf(X,T,rhosol, 'edgecolor','none');
% axis([-L/2 L/2 0 tmax]);view(90,-90);grid off   
% tits = ['A=', num2str(A), ', D=',num2str(D),', \rho=',num2str(r0)];
% title(tits)
% %caxis([0,1]);
% colorbar;view(90,-90)
% xlabel x, ylabel Time
% set(gca,'Fontsize',22,'Fontname','Times New Roman')
% set(0,'defaultaxesLinewidth',1.5);
% set(0,'defaultaxesfontsize',22);
% print(gcf,'-dtiff',filename1)
% 
% 
% % figure(3)
% % p1 = plot(tv,mean_rho);
% % set(p1,'Color','blue','LineStyle','-','LineWidth',2)
% % hold on
% % p2 = plot(tv,mean_moment);
% % set(p2,'Color','red','LineStyle','-','LineWidth',2)
% % p3 = plot(tv,mean_square);
% % set(p3,'Color',[0 0.5 0],'LineStyle','-','LineWidth',2)
% 
% 
% figure(2)
%  M = length(tv);
%  filename2 = ['CH_July22_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_T',num2str(tmax),'.gif'];
%   tsol = tv;
% for i = 1:M
%     if (mod(i,1/dt/10)==1) 
%         h1 = plot(x,rhosol(i,:));
%         xlim([-L/2 L/2]);
%         ylim([0,1]);
%          set(h1,'Color','blue','LineStyle','-','LineWidth',2)
%         title_str=['Time = ', num2str(tsol(i))];
%         title(title_str,'position',[0.2,1])
% % text = ['Time = ', num2str(tsol(i))];
% % annotation('textbox',[0.2 0.1 0.5 0.5],'string',text,'fontsize',20,'color','m','edgecolor','none')
%        xlabel 'x'
%        ylabel '\rho'
%         set(gca,'Fontsize',14,'Fontname','Times New Roman')
%        drawnow;
%       im = frame2im(getframe(gcf));  % Get current figure
%       [AA, map] = rgb2ind(im,256);
%            if (i==1)
%                imwrite(AA,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.0001);
%            else
%                imwrite(AA,map,filename2,'gif','WriteMode','append','DelayTime',0.0001);
%            end
%     end
% end
% 
% 
%         figure(3);clf
%         usol(:,N+1) = usol(:,1);
%         filename3 = ['CH_July22_rhou_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_T',num2str(tmax),'.tif'];
%         
%         yyaxis left
%         h1 = plot(x,rhosol(end,:));
%         set(h1,'Color','blue','LineStyle','-','LineWidth',2);
%         xlabel('x')
%          ylabel('Density \rho_0')
%         
%         yyaxis right
%         h2= plot(x,usol(end,:));
%         set(h2,'Color','red','LineStyle','-','LineWidth',2);
%         ylabel('Displacement u')
%         ymax = max(usol(end,:))
%          ymin = min(usol(end,:))
%         ylim([ymin ymax])
%         title_str=['Time = ', num2str(tsol(end))];
%         title(title_str,'position',[0.2,1])
%         set(0,'defaultlineLinewidth',1.5);
%         set(0,'defaultaxesLinewidth',1.5);
%         set(gca,'Fontsize',18,'Fontname','Times New Roman')
%         print(gcf,'-dtiff',filename3)
%        

%         figure(4)
%         filename4 = ['CHmech1Du_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'.tif'];
%         
%         h2 = plot(x,usol(end,:));
%         xlim([-L/2 L/2]);
%         set(h2,'Color','red','LineStyle','-','LineWidth',2)
%         title(title_str,'position',[0.2,1])
%           set(0,'defaultlineLinewidth',1.5);
%           set(0,'defaultaxesLinewidth',1.5);
%          set(gca,'Fontsize',18,'Fontname','Times New Roman')
%          xlabel 'X'
%          ylabel 'Displacement u'
%          print(gcf,'-dtiff',filename4)
     
%         figure(5)
%         filename5 = ['CHmech1Dv_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'.tif'];
%         vsol(:,N+1) = vsol(:,1);
%         h2 = plot(x,vsol(end,:));
%         xlim([-L/2 L/2]);    
%         set(h2,'Color','red','LineStyle','-','LineWidth',2)
%         title(title_str,'position',[0.2,1])
%         set(0,'defaultlineLinewidth',1.5);
%         set(0,'defaultaxesLinewidth',1.5);
%         set(gca,'Fontsize',18,'Fontname','Times New Roman')
%         xlabel 'X'
%         ylabel 'Velocity v'
%         print(gcf,'-dtiff',filename5)

 
 










