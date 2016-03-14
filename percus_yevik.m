function [plot_r, rdf] = percus_yevik( v_f )
%This function is used to calculate the standard percus_yevik solution, the
%solution only valid in the interval (d,2d), the input value is volume
%fraction
        p1=(1-v_f)*(1-v_f);
        p2=6*v_f*(1-v_f);
        p3=18*v_f*v_f;
        p4=-12*v_f*(1+2*v_f);
        p=[p1,p2,p3,p4];
        t = roots(p);

        times_r = 21;
        delta_r = 0.05;
        plot_r = linspace(1,2,times_r);

        rdf=zeros(times_r,1);
        for i = 1:times_r
            r = 1.0 + (i-1)*delta_r;
            rdf(i) = g_1(r,t,v_f);
        end
end
function g1 = g_1(r,t,vf)
g1=0.0;
for ii=1:3
    g1=g1+t(ii)*Lt(vf,t(ii))*exp(t(ii)*(r-1))/St(vf,t(ii));
end
g1 = g1/(12*vf*r);
end
function [L_t] = Lt(v_f,t)
L_t = 12*v_f*((1+0.5*v_f)*t + 1+2*v_f);
end

function [S_t_prime] = St(v_f, t)
S_t_prime = 3*(1-v_f)*(1-v_f)*t*t + 12*v_f*(1-v_f)*t +18*v_f*v_f;
end
