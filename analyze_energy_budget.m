function[time, norm_t, norm_f, Fhx, Fhz, Fix, Fiz, Fx, Fz, PPE, PKE, FKE, Ex, Ez, PPE_decay, PKE_decay, FKE_decay, epsilon_fluid, epsilon_part] = calculate_energy_budget(casename, savename, time, a, g, rho_p, theta_p);
%clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the potential energy of particles, kinetic energy of particles and fluid
%Calculate the decay rate of each of them,using one-side forward euler
%Calculate the dissipation due to particles and fluid;
%       for particles: is calculated by multiply velocity and hydrodynamic force
%       for fluid: is calculated by 2\mu SijSij; velocity of particle phase is 
%                  represented by solid body 
%Notice:        Particles will not be all on the ground at last time step, Notice
%               the normalization of energy.
%               The normalization of decay rate and dissipation is \Delta E/norm_t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%casename = '';
[~, xe, ~, dx, ~, ~, ~, dy, zs, ~, ~, dz] = calculate_domain_info(casename);
%information regarding the simulations
%a = 1; g = 10;
rho_f = 1;nu = 1;
%rho_p = 2.6; theta_p = 0;
%t1 = 0:0.1:2.5; t2 = 2.8:0.1:2.9;
%time = [t1, t2];
%savename = '_as=1';
tmp = ['force',savename,'.mat'];
savedata{1} = tmp;
tmp = ['energy_budget',savename,'.mat'];
savedata{2} = tmp;

v_p = 4*pi*a^3/3; m_f = dx*dy*dz*rho_f;

%calculate initial potential energy
[x, y, z] = cgns_read_part_position(casename,0);
Np = length(x);
E_ref = (rho_p - rho_f)*v_p*a*g*Np;
E0 = sum((xe - x)*sin(theta_p) + (z - zs)*cos(theta_p)) * (rho_p - rho_f) * v_p * g - E_ref;
norm_t = sqrt(2*a/(g*cos(theta_p)));
norm_f = (rho_p - rho_f)*v_p*g*cos(theta_p);

%accumulate the statistics
for tt = 1:Nt
        fprintf('time is %f\n',time(tt));
        [x, y, z] = cgns_read_part_position(casename,time(tt));
        [u, v, w] = cgns_read_part_vel(casename,time(tt));

        [Fhx, Fhy, Fhz] = cgns_read_part_force_hydro(casename, time(tt));
        Hfx(tt,1) = mean(Fhx);
        Hfz(tt,1) = mean(Fhz);
        [Fix, Fiy, Fiz] = cgns_read_part_force_interaction(casename,time(tt));
        Ifx(tt,1) = mean(Fix);
        Ifz(tt,1) = mean(Fiz);
        [F1, F2, F3] = cgns_read_part_force_total(casename,time(tt));
        Fx(tt,1) = mean(F1);
        Fz(tt,1) = mean(F3);

        [U, V, W] = cgns_read_flow_vel(casename, time(tt));
        [phase] = cgns_read_flow_phase(casename,time(tt));
        phase(phase~=-1) = 0;
        phase(phase==-1) = 1;
        phase = double(phase);
        [FUY, FUX, FUZ] = gradient(U,dy,dx,dz);
        [FVY, FVX, FVZ] = gradient(V,dy,dx,dz);
        [FWY, FWX, FWZ] = gradient(W,dy,dx,dz);

        %%dissipation of fluid is calculated as the velocity of particle is the solid body velocity
        epsilon_fluid(tt,1) = nu * sum(sum(sum(phase.*(2*FUX.*FUX  + 2*FVY.*FVY + 2*FWZ.*FWZ + FUY.*FUY + FVX.*FVX + 2*FUY.*FVX + ...
                              FUZ.*FUZ + FWX.*FWX + 2*FUZ.*FWX + FVZ.*FVZ + FWY.*FWY + 2*FVZ.*FWY))))*dx*dy*dz;
        epsilon_part(tt,1) = -sum(Fhx.*u + Fhy.*v + Fhz.*w);

        %%potential eneregy of particle is calculated as sum(instantenous height),subtract the intial value
        PPE(tt,1) = sum((xe - x)*sin(theta_p) + (z - zs)*cos(theta_p))*(rho_p - rho_f)*v_p*g - E_ref;
        PKE(tt,1) = 0.5 * sum(u.*u + v.*v + w.*w) * v_p * rho_p;
        Ex(tt,1) = 0.5 * sum(u.*u)*v_p * rho_p;
        Ez(tt,1) = 0.5 * sum(w.*w)*v_p * rho_p;

        U = U.*phase; V = V.*phase; W = W.*phase;
        FKE(tt,1) = 0.5 * sum(sum(sum((U.*U + V.*V + W.*W))))*m_f;
        end
end

Fhx = Hfx; Fhz = Hfz;
Fix = Ifx; Fiz = Ifz;

eval(sprintf('time%s = time;',savename));eval(sprintf('norm_t%s = norm_t;',savename)); eval(sprintf('norm_f%s = norm_f;',savename));
eval(sprintf('Fhx%s = Fhx;',savename));eval(sprintf('Fhz%s = Fhz;',savename));
eval(sprintf('Fix%s = Fix;',savename));eval(sprintf('Fiz%s = Fiz;',savename));
eval(sprintf('Fx%s = Fx;',savename));eval(sprintf('Fz%s = Fz;',savename));
save(savedata{1},'time_*','norm_t_*','norm_f_*','Fhx_*','Fhz_*','Fix_*','Fiz_*','Fx_*','Fz_*');

PPE(:,1) = PPE(:,1)/E0;
PKE(:,1) = PKE(:,1)/E0;
FKE(:,1) = FKE(:,1)/E0;
Ex(:,1) = Ex(:,1)/E0;
Ez(:,1) = Ez(:,1)/E0;
epsilon_fluid(:,1) = epsilon_fluid(:,1)/E0(1,1);
epsilon_part(:,1) = epsilon_part(:,1)/E0(1,1);


for tt = 1:Nt-1
        PPE_decay(tt,1) = (PPE(tt+1,1) - PPE(tt,1))/(time(1,tt+1) - time(1,tt));
        PKE_decay(tt,1) = (PKE(tt+1,1) - PKE(tt,1))/(time(1,tt+1) - time(1,tt));
        FKE_decay(tt,1) = (FKE(tt+1,1) - FKE(tt,1))/(time(1,tt+1) - time(1,tt));
end
eval(sprintf('PPE%s = PPE;',savename)); eval(sprintf('PKE%s = PKE;',savename));eval(sprintf('FKE%s = FKE;',savename));
eval(sprintf('Ex%s = Ex;',savename));eval(sprintf('Ez%s = Ez;',savename));
eval(sprintf('PPE_decay%s = PPE_decay;',savename)); eval(sprintf('PKE_decay%s = PKE_decay;',savename)); eval(sprintf('FKE_decay%s = FKE_decay;',savename));
eval(sprintf('epsilon_fluid%s = epsilon_fluid;',savename)); eval(sprintf('epsilon_part%s = epsilon_part;',savename));
save(savedata{2},'time_*','norm_t_*','PPE_*','PKE_*','Ex_*','Ez_*','FKE_*','PPE_decay_*','PKE_decay_*','FKE_decay_*','epsilon_fluid_*','epsilon_part_*');
end
