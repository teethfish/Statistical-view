function [ s, Cqq ] = calculate_correlation_for_1d_array(qx, qy)
qx = qx - mean(qx);
qy = qy - mean(qy);
Nt = length(qx);
%Nt - 2 is the total length, use Nt/2 to avoid sampling error
s = 0:1:fix(Nt/2);
Ns = length(s);
Cqq = zeros(Ns,1);
for ss = 1:Ns
        ds = s(ss);
        pair = Nt - ds;
        vector_i = qx(1:pair);
        vector_j = qy(1+ds:pair+ds);
        Cqq(ss,1) = sum(vector_i.*vector_j)/pair;
        %for pp = 1:pair
        %        tmp_1 = pp;
        %        tmp_2 = tmp_1 + ds;
        %        Cqq(ss,1) = Cqq(ss,1) + qx(tmp_1) * qy(tmp_2);
        %end
        %Cqq(ss,1) = Cqq(ss,1)/pair;
end
Cqq = Cqq/Cqq(1,1);
end
