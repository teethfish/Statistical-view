function [nq] = dam_breaking_calculate_part_2d_fft(q, x, z, Hx, Hz, M, L)
% v is the quanlity to expansion(velocity,number density...)
% x, z is the two direction for expansion;
% Hx, and Hz is the half length in two direction;
% M, L is the expansion coefficents in x,z direction.

nq = zeros(2*M+1, 2*L+1);
km = (-M:M)'*pi/Hx;
kl = (-L:L)'*pi/Hz;
V = 2*Hx*2*Hz;
for m = -M:M
    for l = -L:L
        nq(m+M+1, l+L+1) = sum(exp(1i*(km(m+M+1)*x + kl(l+L+1)*z)).*q)/V;
    end
end

end
