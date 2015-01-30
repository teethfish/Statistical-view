function pn = legendre(n,z)
    pn=0;
    if n==0
        pn=1;
    end
    if n==1
        pn=z;
    end
    if n==2
        pn=0.5*(3*z*z-1);
    end
    if n==3
        pn=0.5*(5*z*z*z-3*z);
    end
    if n==4
        pn=0.125*(35*z*z*z*z-30*z*z+3);
    end
    if n==5
        pn=0.125*(63*z*z*z*z*z-70*z*z*z+15*z);
    end
    if n==6
        pn=(1/16)*(231*z*z*z*z*z*z-315*z*z*z*z+105*z*z-5);
    end
end