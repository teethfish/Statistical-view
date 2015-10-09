clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the particle potential energy, when the ratio decreases to a certian value, 
%record the time and crecord the wave front position and velocity
%(furthest particle and velocity of that furthest particle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

casename{1} = '';
casename{2} = '';
Ncase = length{casename};

%when PPE decreases to critial_PPE% of initial PPE
critical_PPE = 0.2; 
%how many particles to track
Nm = 30;

[xs, xe, ~, dx, ys, ~, ~, dy, zs, ~, ~, dz] = calculate_domain_info(casename{1});
%set the time
t = 0:0.1:30;
Nt = length(t);

%set the simulation property
a = 1; g = 10;
rho_p{1} = 5; rho_f{1} = 1; theta_p{1} = 0;

v_p = 4*pi*a^3/3;  m_f = dx*dy*dz*rho_f;
for i = 1:Ncase
  m_p{i} = v_p*rho_p{i};
  norm_t{i} = sqrt(2*a/(g*cos(theta_p{i})));
end  

for i = 1:Ncase 
  [x, y, z] = cgns_read_part_position(casename{1},0);
  Np = length(x);

  %reference potential energy, when all the particles are on the ground
  E_ref{i} = (rho_p{i} - rho_f)*v_p*a*g*Np;

  %initial potential energy
  E0{i} = sum((xe - x)*sin(theta_p{i}) + (z - zs)*cos(theta_p{i})) * (rho_p{i} - rho_f) * v_p * g - E_ref{i};
end

for i = 1:Ncase
        [~, ~, z0] =cgns_read_part_position(casename{i},0);
        Np = length(z0);
        [~,tmp] = sort(z0,'descend');
        I = tmp(1:Nm);
        I = sort(I,'ascend');
        for j = 1:Np
                for k = 1:Nm
                        if j==I(k)
                                index(j,i) = 1;
                                break;
                        else
                                index(j,i) = 0;
                        end
                end
        end
end

%record the wave from position and wave front velocity by the time when particle potential decrease to critical value
for tt = 1:Nt
        fprintf('time is %f\n',t(tt));
        for i = 1:Ncase
          [x, y, z] = cgns_read_part_position(casename{i},t(tt));
          [u, v, w] = cgns_read_part_vel(casename{i},t(tt));
          PPE(tt,i) = (sum((xe - x)*sin(theta_p{i}) + (z - zs)*cos(theta_p{i}))*(rho_p{i} - rho_f)*v_p*g - E_ref{i})/E0{i};
          [tmp, index] = max(x);
          wave_front(tt,i) = tmp - xs;
          wave_front_speed(tt,i) = max(u);
          avg_vel(tt,i) = sum(u.*index(:,i))/Nm;
          if PPE(tt,1) <= critical_PPE
              te = t(tt);
              break;
          end
        end  
end

time = 0:0.1:te;

%rename variables and save in .mat file
time_rho_5 = time;
norm_t_rho_5 = norm_t;
wave_front_rho_5 = wave_front;
wave_front_speed_rho_5 = wave_front_speed;
save('wave_front_rho_5.mat','time_rho_5','norm_t_rho_5','wave_front_rho_5','wave_front_speed_rho_5');

