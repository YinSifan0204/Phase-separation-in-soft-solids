%%%%%%%%%%%%%%%%%%%
%  Analysis of 1D simulation        %
%       gamma = 0.01                    %
%                                                  %
%        Sifan Yin                            %
%      04/22/2020  
%     07/22/2020
%%%%%%%%%%%%%%%%%%%

% newcolorbar=colormap;
% save mycolor newcolorbar;

clc;clear all;close all;

A=5
r0 = 0.2
D= 0.01;
gamma = 0.1;
%paraname = ['Parameters1D0420_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'.mat'];
paraname = ['Parameters1D_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'.mat']
load(paraname)
gamma = para(4); 
dt = para(5);
tmax = para(6);
N = para(7);
L = 10;

load('BlueGrayColormap.mat');

%dataname =  ['Data1D0420_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat'];
dataname =  ['Data1D_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat'];
load(dataname);
M = length(tv);

% Extract the fft 
dx = L/N;
x = L/N*[-N/2:N/2-1];
kX = (2*pi)/L*[0:N/2-1 -N/2:-1];
rho_fft = real(fft(rhosol(M,:)));
[maxrho,index]=max(rho_fft);
rho_fft(index)=0;
[maxrho,index]=max(rho_fft)
kmax = kX(index)
lambda=2*pi/kmax




% for i = 1:M
%         mean_rho(i) = mean(rhosol(i,:));
%         mean_moment(i) = mean(rhosol(i,:).*x);
%         mean_square(i)= mean(rhosol(i,:).*((x - mean_moment(i)).^2));
% end
figure('color',[1 1 1]); % White backgraoud
set(gcf, 'position',[200 100 1200 1200]);


filename1 = ['CH1D_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),...
                      '_T',num2str(tmax),'_rho.tif'];                            
 x(N+1) = L/2; 
 %tv = tv(1:5001)
 [X,T] = meshgrid(x,tv);
 rhosol(:,N+1) = rhosol(:,1);
 %rhosol = rhosol(1:5001,:);
surf(X,T,rhosol, 'edgecolor','none');
grid on
%axis([-L/2 L/2 0 tmax/2]);view(90,-90);grid off   
axis([-L/2 L/2 0 tmax]);view(90,-90);grid off   
%axis equal
%tits = ['A=', num2str(A), ', D=',num2str(D),', \rho=',num2str(r0)];
%tits = 'Cell density \rho_0'
%title(tits)
%caxis([0,1]);
colormap(BlueGrayColormap);
% colorbar;
% print(gcf,'-dtiff','BlueGrayColorbar')

view(90,-90)
axis off
%xlabel x, ylabel Time
% set(gca,'Fontsize',18,'Fontname','Times New Roman')
% set(gca,'Fontsize',22,'Fontname','Times New Roman')
% set(0,'defaultaxesLinewidth',1.5);
% set(0,'defaultaxesfontsize',22);
print(gcf,'-dtiff',filename1)

% figure(2)
% filename2 = ['CH1D0420_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),...
%                       '_T',num2str(tmax),'_u.tif'];
%  usol(:,N+1) = usol(:,1);
% surf(X,T,usol, 'edgecolor','none');
% grid on
% axis([-L/2 L/2 0 tmax]);view(90,-90);grid off   
% %tits = ['A=', num2str(A), ', D=',num2str(D),', \rho=',num2str(r0)];
% tits = 'Displacement u'
% title(tits)
% %caxis([0,1]);
% colorbar;view(90,-90)
% xlabel x, ylabel Time
% set(gca,'Fontsize',18,'Fontname','Times New Roman')
% set(gca,'Fontsize',22,'Fontname','Times New Roman')
% set(0,'defaultaxesLinewidth',1.5);
% set(0,'defaultaxesfontsize',22);
% print(gcf,'-dtiff',filename2)
% 
% figure(2)
% p1 = plot(tv,mean_rho);
% set(p1,'Color','blue','LineStyle','-','LineWidth',2)
% hold on
% p2 = plot(tv,mean_moment);
% set(p2,'Color','red','LineStyle','-','LineWidth',2)
% p3 = plot(tv,mean_square);
% set(p3,'Color',[0 0.5 0],'LineStyle','-','LineWidth',2)
% filename3 = ['CH_Apr14ps_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'.tif']; 
% xlabel 'Time'
% ylabel 'P1,P2,P3'
% set(gca,'Fontsize',14,'Fontname','Times New Roman')
% legend('P1','P2','P3')
% print(gcf,'-dtiff',filename3)


% figure(3)
%  filename3 = ['CH_July29_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_T',num2str(tmax),'_rho.gif'];
%   tsol = tv;
% for i = 1:M
%     if (mod(i,0.1/dt)==1) % Every 0.1s outputs an image 
%         h1 = plot(x,rhosol(i,:));
%         xlim([-L/2 L/2]);
%         ylim([0,1]);
%          set(h1,'Color','blue','LineStyle','-','LineWidth',2)
%         title_str=['Time = ', num2str(tsol(i))];
%         title(title_str,'position',[0.2,1.2])
% % text = ['Time = ', num2str(tsol(i))];
% % annotation('textbox',[0.2 0.1 0.5 0.5],'string',text,'fontsize',20,'color','m','edgecolor','none')
%        xlabel 'x'
%        ylabel '\rho'
%         set(gca,'Fontsize',14,'Fontname','Times New Roman')
%        drawnow;
%       im = frame2im(getframe(gcf));  % Get current figure
%       [AA, map] = rgb2ind(im,256);
%            if (i==1)
%                imwrite(AA,map,filename3,'gif','LoopCount',Inf,'DelayTime',0.01);
%            else
%                imwrite(AA,map,filename3,'gif','WriteMode','append','DelayTime',0.01);
%            end
%     end
% end



% figure(4)
%  filename4 = ['CH1D0420_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_T',num2str(tmax),'_rhou.gif'];
%   tsol = tv;
% for i = 1:M
%  if (mod(i,0.01/dt)==1) % Every 0.1s outputs an image 
%         h1 = plot(x,rhosol(i,:));
%         set(h1,'Color','blue','LineStyle','-','LineWidth',2);
%         xlim([-L/2 L/2]);
%         ylim([-1,1.1]);
%          set(h1,'Color','blue','LineStyle','-','LineWidth',2)
%          title_str=['Time = ', num2str(tsol(i))];
%          title(title_str,'position',[0.2,1.1])
% % text = ['Time = ', num2str(tsol(i))];
% % annotation('textbox',[0.2 0.1 0.5 0.5],'string',text,'fontsize',20,'color','m','edgecolor','none')
%         xlabel 'x'
%         ylabel('Cell density \rho')
% 
%         yyaxis right
%         h2= plot(x,usol(i,:));
%         set(h2,'Color','red','LineStyle','-','LineWidth',2);
%         ylabel('Displacement u')
%         ylim([-1 3])
% %         ymax = max(usol(i,:))
% %          ymin = min(usol(i,:))
%         %ylim([ymin ymax])
%         set(gca,'Fontsize',22,'Fontname','Times New Roman')
%         
%         
%        drawnow;
%       im = frame2im(getframe(gcf));  % Get current figure
%       [AA, map] = rgb2ind(im,256);
%            if (i==1)
%                imwrite(AA,map,filename4,'gif','LoopCount',Inf,'DelayTime',0.01);
%            else
%                imwrite(AA,map,filename4,'gif','WriteMode','append','DelayTime',0.01);
%            end
%            clf;
%     end
% end






%         figure(5);clf
%         toutput = tmax;
%         index = round(length(tsol))
%         filename5 = ['CH1D0420_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'_T',num2str(toutput),'_ru.tif'];
%         yyaxis left
%         
%          h1 = plot(x,rhosol(index,:));
%          set(h1,'Color','blue','LineStyle','-','LineWidth',2);
%          xlabel('x')
%          ylabel('Density \rho')
%          ylim([-1,1.1])
%          
%         yyaxis right
%         h2= plot(x,usol(index,:));
%         set(h2,'Color','red','LineStyle','-','LineWidth',2);
%         ylabel('Displacement u')
%         ymax = max(usol(index,:))
%          ymin = min(usol(index,:))
%         %ylim([ymin ymax])
%         ylim([-1 3])
%         title_str=['Time = ', num2str(tsol(index))];
%         title(title_str,'position',[0.2,1.1])
%         set(0,'defaultlineLinewidth',1.5);
%         set(0,'defaultaxesLinewidth',1.5);
%         set(gca,'Fontsize',22,'Fontname','Times New Roman')
%         print(gcf,'-dtiff',filename5)
        
        
        
        
        
%         figure(6);clf
%         toutput = tmax;
%         index = round(length(tsol))
%         rhosol_temp = rhosol(index,:);
%         rhoh = fft(rhosol_temp)/N;
%          kx = (2*pi/L)*[0:N/2 -N/2:-1];
%          h1 = plot(kx,rhoh);
%          set(h1,'Color','blue','LineStyle','-','LineWidth',2);
%          xlabel('x')
%          ylabel('Density \rho')
%          ylim([-1,1.1])
%          
%         yyaxis right
%         h2= plot(x,usol(index,:));
%         set(h2,'Color','red','LineStyle','-','LineWidth',2);
%         ylabel('Displacement u')
%         ymax = max(usol(index,:))
%          ymin = min(usol(index,:))
%         %ylim([ymin ymax])
%         ylim([-1 3])
%         title_str=['Time = ', num2str(tsol(index))];
%         title(title_str,'position',[0.2,1.1])
%         set(0,'defaultlineLinewidth',1.5);
%         set(0,'defaultaxesLinewidth',1.5);
%         set(gca,'Fontsize',22,'Fontname','Times New Roman')
%        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
%        

%         figure(4)
%         filename4 = ['CH_July23_L10_N',num2str(N),'_A', num2str(A), '_D',num2str(D),'_rho',num2str(r0),'.tif'];
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










