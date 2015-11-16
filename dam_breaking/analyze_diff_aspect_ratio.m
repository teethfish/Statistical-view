casename{1} = 
time{1} = 
savename{1} = 
casename{2} = 
time{2} = 
savename{2} = 
casename{3} = 
time{3} =
savename{3} = 
Ncase = length{casename};

a = 1; 
g = 10;
rho_p = 3;
theta_p = 0;
critical_PPE = 0.2;
Nm = 40;
for i = 1:Ncase
  [~, norm_t{i}, norm_f{i}, Fhx{i}, Fhz{i}, Fix{i}, Fiz{i}, Fx{i}, Fz{i}, PPE{i}, PKE{i}, FKE{i}, ...
  Ex{i}, Ez{i}, PPE_decay{i}, PKE_decay{i}, FKE_decay{i}, epsilon_fluid{i}, epsilon_part{i}] = calculate_energy_budget(casename{i}, savename{i}, time{i}, a, g, rho_p, theta_p);
  [~, u_mean{i}, w_mean{i},x_mean{i}, z_mean{i}] = calculate_part_distri_layer(casename{i}, time{i}, savename{i});
  [~,~,wave_front{i},wave_front_speed{i},avg_vel{i}] = calculate_wave_front(casename{i}, time{i}, savename{i}, rho_p, theta_p, a, g, critical_PPE, Nm);
end

save('force.mat',time, norm_t, norm_f, Fhx, Fhz, Fix, Fiz, Fx, Fz);
save('energy_budget.mat', time, norm_t, PPE, PKE, FKE, Ex, Ez, PPE_decay, PKE_decay, FKE_decay, epsilon_fluid, epsilon_part);
save('part_distri.mat',time, norm_t, u_mean, w_mean, x_mean, z_mean);
save('wave_front.mat',time,norm_t, wave_front, wave_front_speed, avg_vel);




