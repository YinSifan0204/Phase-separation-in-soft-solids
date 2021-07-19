%%%%%%%%%%%%%%%%%%%%%%%%
%    Parameters list
%     Sifan Yin 
%    8/28/2020
%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;
A= 10;
r0 = 0.5;

D= 0.01;
gamma = 0.01;

dt = 0.05;
ist = 20;  % Output interval
tmax = 100;
N = 512;
 L = 10;
para2D = [A, D, r0, gamma,dt,tmax,N,ist,L];
paraname = ['Parameters2D_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'_gamma',num2str(gamma),'.mat']
save(paraname,'para2D')