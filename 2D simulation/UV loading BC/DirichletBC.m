function [dleft,dright] = DirichletBC(N,L)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
 n = 0:N;
 o = ones(1,N+1);
 dleft = (-o).^n;  % dleft*uch=f(x0)
 dright = o ;  %dright*uch = f(xN)
end

