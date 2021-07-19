%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Cahn-Hilliard Equ with mechanical effect          %
%               Linear stability analysis                               %
%       Dispersion Relation & Phase diagram                %
%                                                                                   %  
%                           Sifan Yin                                          %
%                                        diff = 0
%                           08/03/2020                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; fclose('all');

%% 1. Parameters
% Mechanical properties
A    =  5;
rho0 = 0.2;
E_a =   0.1;
zeta_a = 0.1;
KN  = 10;

%% 1. Parameters
% Basic quantities

E_p      = 1;
zeta_p = 1;
eta      = 1;
Diff     = 0.01;%0.02;     %if eta=0.1, then Diff=Diff_bar*eta=0.01;
%title_str = strcat('Diffusivity = ',num2str(Diff));
K = 2;
m = 2;
alpha      =      -1;
beta        =      1;
gamma  =    0.1;
E_a =   0.1;
zeta_a = 0.1;
filename = ['Dispersion relation2021','_A', num2str(A), '_D',num2str(Diff),'_rho',num2str(rho0),'.tif'];

 
% % Dimensionless parameters
alpha      =   alpha/E_p;
beta       =   beta/E_p;
gamma  =    gamma*eta/(E_p*zeta_p);

Diff = Diff*eta;
E_a =   E_a/E_p;
zeta_a = zeta_a/zeta_p;
E_rho0 = E_a*rho0+(1-rho0);
zeta_rho0 = zeta_a*rho0+(1-rho0);



%% Dispersion relation
kn = 0:0.0001:KN;
abc = alpha+3*beta*rho0^2+gamma*kn.^2;
Denominator = 1+zeta_rho0*kn.^2;

fp0  =  m*rho0^(m-1)/(1+K*rho0^m)^2;

bk = kn.^2.*(  Diff*abc + (E_rho0-rho0*fp0*A)./Denominator );
ck = Diff*kn.^4.*( E_rho0*abc - A^2*fp0.^2 )./Denominator;
DELTA = bk.^2-4*ck;
r = -bk./2;
for ii = 1:length(kn)
    if DELTA(ii)>=0
        omega(ii) = 0;
        w = sqrt(DELTA(ii))/2;
        lambda1(ii) = r(ii) + w;
        lambda2(ii) = r(ii)  - w;
    else
        omega(ii) = sqrt(-DELTA(ii))/2;
        lambda1(ii) = r(ii);
        lambda2(ii) = r(ii);
    end
end
 y0 = zeros(1,length(kn));
% figure(1)
% plot(kn,DELTA,'r-')
figure(2)
%title_str= ['A=', num2str(A), ', D=',num2str(Diff),', \rho_0=',num2str(rho0)];
plot(kn,y0,'k--',kn,omega,'r--',kn,-omega,'r--',kn,lambda1,'b-',kn,lambda2,'b-','LineWidth', 2)
%legend('zero','\omega(k_n)','-\omega(k_n)','\lambda_1','\lambda_2','Location','best')
xlim([0 KN])
ylim([-10 2])
set(0,'defaultlineLinewidth',1.5);
set(0,'defaultaxesLinewidth',1.5);
set(gca,'Fontsize',18,'Fontname','Times New Roman')
xlabel k_n
ylabel \lambda
%title(title_str)
set(gca,'Fontsize',18,'Fontname','Times New Roman')
print(filename,'-dtiff')

