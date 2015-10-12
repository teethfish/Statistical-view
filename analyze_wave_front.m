function[time,norm_t,wave_front,wave_front_speed,avg_vel] = calculate_wave_front(casename, t, savename, rho_p, theta_p, a, g, critical_PPE, Nm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the particle potential energy, when the ratio decreases to a certian value, 
%record the time and crecord the wave front position and velocity
%(furthest particle and velocity of that furthest particle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%casename = '';
%when PPE decreases to critial_PPE% of initial PPE
%critical_PPE = 0.2; 
%how many particles to track
%Nm = 30;

[xs, xe, ~, dx, ys, ~, ~, dy, zs, ~, ~, dz] = calculate_domain_info(casename);
savedata = ['wave_front',savename,'.mat'];
%t = 0:0.1:30;
dt = t(2) - t(1);
Nt = length(t);

%set the simulation property
rho_f = 1;
%a = 1; g = 10;
%rho_p = 5;  theta_p= 0;

v_p = 4*pi*a^3/3;  m_f = dx*dy*dz*rho_f;
m_p = v_p*rho_p;
norm_t = sqrt(2*a/(g*cos(theta_p)));

[x, y, z] = cgns_read_part_position(casename,0);
Np = length(x);
%reference potential energy, when all the particles are on the ground
E_ref = (rho_p - rho_f)*v_p*a*g*Np;
%initial potential energy
E0 = sum((xe - x)*sin(theta_p) + (z - zs)*cos(theta_p)) * (rho_p - rho_f) * v_p * g - E_ref;

%Find initially highest Nm particles and set an array index for it
[~, ~, z0] =cgns_read_part_position(casename,0);
Np = length(z0);
[~,tmp] = sort(z0,'descend');
I = tmp(1:Nm);
I = sort(I,'ascend');
for j = 1:Np
    for k = 1:Nm
       if j==I(k)
          index(j,1) = 1;
          break;
        else
          index(j,1) = 0;
        end
    end
end

%record the wave from position and wave front velocity by the time when particle potential decrease to critical value
for tt = 1:Nt
        fprintf('time is %f\n',t(tt));
        [x, y, z] = cgns_read_part_position(casename,t(tt));
        [u, v, w] = cgns_read_part_vel(casename,t(tt));
        PPE(tt,1) = (sum((xe - x)*sin(theta_p) + (z - zs)*cos(theta_p))*(rho_p - rho_f)*v_p*g - E_ref)/E0;
        [tmp, index] = max(x);
        wave_front(tt,1) = tmp - xs;
        wave_front_speed(tt,1) = max(u);
        avg_vel(tt,1) = sum(u.*index(:,1))/Nm;
        if PPE(tt,1) <= critical_PPE
          te = t(tt);
          break;
        end
    end  
end

time = t(1):dt:te;
%rename variables and save in .mat file
eval(sprintf('time%s = time;',savename));
eval(sprintf('norm_t%s = time;',savename));
eval(sprintf('wave_front%s = wave_front;',savename));
eval(sprintf('wave_front_speed%s = wave_front_speed;',savename));
eval(sprintf('avg_vel%s = avg_vel;',savename));
save(savedata,'time_*','norm_t_*','wave_front_*','wave_front_speed_*','avg_vel_*');

end
