function [dleft,dright] = DirichletBC(N,L)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
 n = 0:N;
 o = ones(1,N+1);
 dleft = (-o).^n;  % dleft*uch=f(x0)
 dright = o ;  %dright*uch = f(xN)
end

