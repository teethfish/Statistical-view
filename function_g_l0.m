function [gl_x, gl_y, gl_z ] = function_g_l0(NP, l, d, dgr, x, y, z, mg, Hx, Hy, Hz, xs, ys, zs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function is used to get the gl0 of the radial distribution function
%NP:    the total number of particles
%l:     order of gl0
%d:     diameter of the particle
%dgr:   the interval of r
%x,y,z: the coordinate of particle in three directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = x;
    p = y;
    n = z;  
    gl_x = zeros(mg,1);
    gl_y = zeros(mg,1);
    gl_z = zeros(mg,1);
    V = Hx*Hy*Hz;
    for i = 1:NP
        for j =1:NP
            if x(i) < (xs+Hx/2)
                if x(j)-x(i) > Hx/2
                    m(j)=x(j)-Hx;
                else
                    m(j)=x(j);
                end
            end
            if x(i) > (xs+Hx/2)
                if x(i)-x(j) > Hx/2
                    m(j)=x(j)+Hx;
                else
                    m(j)=x(j);
                end
            end
            if y(i) < (ys+Hy/2)
                if y(j)-y(i) > Hy/2
                    p(j)=y(j)-Hy;
                else
                    p(j)=y(j);
                end
            end
            if y(i) > (ys+Hy/2)
                if y(i)-y(j) > Hy/2
                    p(j)=y(j)+Hy;
                else
                    p(j)=y(j);
                end
            end
            if z(i) < (zs+Hz/2)
                if z(j)-z(i) > Hz/2
                    n(j)=z(j)-Hz;
                else
                    n(j)=z(j);
                end
            end
            if z(i) > (zs+Hz/2)
                if z(i)-z(j) > Hz/2
                    n(j)=z(j)+Hz;
                else
                    n(j)=z(j);
                end    
            end  
            d_two = sqrt((m(j)-x(i))*(m(j)-x(i))+(p(j)-y(i))*(p(j)-y(i))+(n(j)-z(i))*(n(j)-z(i))); %Real distance between two particles
            tmp_d = fix((d_two/d - 1)/dgr)+1; %Change the distance into an integer
            if tmp_d > 0 && tmp_d <= mg
                cos_theta_x = (m(j)-x(i))/d_two;
                cos_theta_y = (p(j)-y(i))/d_two;
                cos_theta_z = (n(j)-z(i))/d_two;
                gl_x(tmp_d) = gl_x(tmp_d) + legendre(l,cos_theta_x)/d_two/d_two;
                gl_y(tmp_d) = gl_y(tmp_d) + legendre(l,cos_theta_y)/d_two/d_two;
                gl_z(tmp_d) = gl_z(tmp_d) + legendre(l,cos_theta_z)/d_two/d_two;
            end
        end
    end
    gl_x = gl_x*V*sqrt((2*l+1)/(4*pi))/(NP*(NP-1)*dgr);
    gl_y = gl_y*V*sqrt((2*l+1)/(4*pi))/(NP*(NP-1)*dgr);
    gl_z = gl_z*V*sqrt((2*l+1)/(4*pi))/(NP*(NP-1)*dgr);
end

