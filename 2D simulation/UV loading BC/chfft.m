function uch= chfft(u)
% Chebyshev expansion truncation coefficients
M = length(u);
N = M-1;

 v = zeros(2*N,1);
 v(1:N+1) = u(end:-1:1); % left-mirror
 v(N+2:end) = u(2:N);
 a = real(fft(v)/N);
 a = a(1:N+1);
 a(1) = a(1)/2;
 a(N+1) = a(N+1)/2;
 uch = a;
end

