%%%%%%%%%%%%%%%%%%%%%%%%
%    Parameters list
%     Sifan Yin 
%    7/23/2020
%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;
A=4
r0 = 0.6
D= 0.01;
gamma = 0.1;

dt = 0.001;
ist = 100;  % Output interval
tmax = 20;
N =512;
para = [A, D, r0, gamma,dt,tmax,N,ist];
paraname = ['Parameters1D0817_D',num2str(D),'_A',num2str(A),'_rho',num2str(r0),'.mat']
save(paraname,'para')