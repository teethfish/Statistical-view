function [nq] = calculate_part_fft(v, y, H, M)
%This function is used to do the scalar expansiton of a particle in 1d
%       H is the half length in the direction
%       M is the order of expansiton of fourier expansion(the total expansion number is 2*M+1)
%       nq is the fft coefficents in complex form
V  = 2*H;
nq = zeros(2*M+1,1);
km = (-M:M)'*pi/H;

for m = -M:M
    nq(m+M+1) = sum(exp(1i*km(m+M+1)*y).*v)/V;
end
%nq_cos = real(nq); nq_sin = imag(nq);
end
