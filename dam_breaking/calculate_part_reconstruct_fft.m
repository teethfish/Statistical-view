function [NQ, beta] = calculate_part_reconstruct_fft(yy, M, nq, H, Np, a, if_cal_beta)
%This function is used to reconstruct the continous field from the fft coefficients
%   yy is coordinates in the expansion directory
%   M is the expansion order
%   nq is the fft coeffients
%   H is the half lengt
%   Np is the particle number
%   a is the particle radius
%   beta is the volume fraction(different from number density)

NUM = length(yy);
km = (-M : M)'*pi/H;
cesaro_coeff = 1-abs((-M : M))'/(M+1);
NQ = zeros(NUM,1);
for j = 1 : NUM
    NQ(j,1) = sum(cesaro_coeff.*nq.*exp(-1i*km*yy(j)));
end
NQ = real(NQ)*2*H/Np;
beta=zeros(NUM,1);

%reconstruct the volume fraction beta, if the input scalar is number density
if(if_cal_beta == 1)
    for j=1:NUM
        for m=-M:M
            k = m*pi/H;
            if m==0;
                volume_coeff=4*pi*a*a*a/3;
            else
                volume_coeff = 4*pi*(sin(abs(k)*a)-abs(k)*a*cos(abs(k)*a))/(abs(k*k*k));
            end
            %beta(j,1)=beta(j,1)+volume_coeff*nq(m+M+1)*exp(-1i*km*yy(j));
            beta(j,1)=beta(j,1)+(1-abs(m)/(M+1))*volume_coeff*nq(m+M+1)*exp(-1i*k*yy(j));
        end
    end
    %to Normalize the volume fraction, divide Xl*Xz, since this is the result for a slice in y-direction
    %beta=real(beta)/50;
    beta = real(beta);
end

end
