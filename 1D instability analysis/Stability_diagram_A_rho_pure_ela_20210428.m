%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Cahn-Hilliard Equ with mechanical effect          %
%               Linear stability analysis                               %
%       Stability Diagram Figure(1-2)                             %
%                           rho ~ A                                           %  
%                           Sifan Yin                                          %
%                        03/09/2021 pure elasticity                 %
%           Non-dimensionalize equations                      %      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; fclose('all');

%% 1. Parameters

E_p      = 1;
zeta_p = 1;
eta      = 1;
Diff = 0.01;
title_str = strcat('{\it D} = ',num2str(Diff));
K = 2;
m = 2;
alpha      =   -1;
beta        =   1;
gamma  =    0.1;%0.02;
E_a =   0.1;
zeta_a = 0.1;

 
% % Dimensionless parameters
alpha      =   alpha/E_p;
beta       =   beta/E_p;
gamma  =    gamma*eta/(E_p*zeta_p);

Diff = Diff*eta;
E_a =   E_a/E_p;
zeta_a = zeta_a/zeta_p;




kb = 1;
kc=4;
Av = 0:0.01:5;
rhov =0:0.001:1;
[A, rho] = meshgrid(Av,rhov);
rhoCfp = m*A.*rho.^m./(1+K*rho.^m).^2 ;
fp = m*rho.^(m-1)./(1+K*rho.^m).^2;
Erho = E_a*rho+1-rho;


bkn2=Diff*kb^2*(alpha+3*beta*rho.^2+gamma*kb^2)+kb^2*(  Erho - rhoCfp  );
ckn2=Diff*kc^4*(  Erho.*(alpha+3*beta*rho.^2+gamma*kc^2) - A.^2.*fp.^2  );
Delkn2 = bkn2.^2-4*ckn2;


figure('color',[1 1 1]); % White backgraoud

bkn2(:,210:end)=NaN;
[M,hbkn2] = contour(A,rho,bkn2,[0 0]);%,[0 2.282 0 1]);
set(hbkn2,'Color','r','LineStyle','-','LineWidth',3)
hold on
ckn2(567:end,:)=NaN;
[M,hckn2] = contour(A,rho,ckn2,[0 0]);%,[0 2.282 0 1]);
set(hckn2,'Color',[0 0.5 0],'LineStyle','-','LineWidth',3)      

set(gca,'Fontsize',22,'Fontname','Times New Roman')
xlabel 'Contractility  \it C'
ylabel 'Cell density \rho_0'
%legend('b=0','c=0','\Delta=0','E-A\rho f(\rho)=0','Location','best')
title(title_str)
set(0,'defaultlineLinewidth',1.5);
set(0,'defaultaxesLinewidth',1.5);
set(0,'defaultaxesfontsize',24);
xlim([0,5])
ylim([0,1])




%% 3/9/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the kn2 corresponding to the largest growth rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kn2m = 1/(2*Diff*gamma)*(   rhoCfp-Erho-Diff*(alpha+3*beta*rho.^2) );
index = find(kn2m<0);
kn2m(index) = NaN;
bkn2m = Diff*kn2m.*(alpha+3*beta*rho.^2+gamma*kn2m)+kn2m.*(  Erho - rhoCfp  );
ckn2m = Diff*kn2m.^2.*( Erho.*(alpha+3*beta*rho.^2+gamma*kn2m) - A.^2.*fp.^2  );
Delm = bkn2m.^2-4*ckn2m;
%Delm(80:end,:)=NaN;
[M,hDelm] = contour(A,rho,Delm,[0 0]);%,[0.4 5 0.6418 1]);
%[M,hDelm] = contour(A,rho,Delm,'ShowText','on');
set(hDelm,'Color','b','LineStyle','-','LineWidth',3)    
 


% end

filename = ['stability_diagram20210428_pureEla_partdel','.tif']
print(gcf,'-dtiff',filename)


% h1 = fimplicit(bk2_old);
% set(h1,'Color','r','LineStyle','-','LineWidth',2,'Marker','o')
% hold on
% h2 = fimplicit(bk2_new);
% set(h2,'Color',[0,0.5,0],'LineStyle','-','LineWidth',2)%,'Marker','o'
% hold on
% h3 = fimplicit(ck2_old);
% set(h3,'Color','b','LineStyle','-','LineWidth',2)
%  h4 = fimplicit(ck2_new);
%  set(h4,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')








