clc; clear all; close all;
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

casename{1} = '';
casename{2} = '';
casename{3} = '';

[~, xe, ~, dx, ~, ~, ~, dy, zs, ~, ~, dz] = calculate_domain_info(casename{1});
a = 1; g = 10;
rho_f = 1;nu = 1;
v_p = 4*pi*a^3/3; m_f = dx*dy*dz*rho_f;

rho_p{1} = 2.6; theta_p{1} = 0;
rho_p{2} = 2.6; theta_p{2} = 0;
rho_p{3} = 2.6; theta_p{3} = 0;
rho_p{4} = 2.6; theta_p{4} = 0;

t1 = 0:0.1:2.5; t2 = 2.8:0.1:2.9; t3 = 3.1:0.1:3.1; t4 = 3.4:0.1:4.7; t5 = 4.9:0.1:5.4; t6 = 5.5:0.1:5.6; t7 = 5.7:0.1:24;
time = [t1, t2, t3, t4, t5, t6, t7];

%calculate initial potential energy
for i = 1:4
        [x, y, z] = cgns_read_part_position(casename{i},0);
        Np = length(x);
        E_ref(i,1) = (rho_p{i} - rho_f)*v_p*a*g*Np;
        E0(i,1) = sum((xe - x)*sin(theta_p{i}) + (z - zs)*cos(theta_p{i})) * (rho_p{i} - rho_f) * v_p * g - E_ref(i,1);
        norm_t{i} = sqrt(2*a/(g*cos(theta_p{i})));
        norm_f{i} = (rho_p{i} - rho_f)*v_p*g*cos(theta_p{i});
end

for tt = 1:Nt
        fprintf('time is %f\n',time(tt));
        for i = 1:4
                [x, y, z] = cgns_read_part_position(casename{i},time(tt));
                [u, v, w] = cgns_read_part_vel(casename{i},time(tt));

                [Fhx, Fhy, Fhz] = cgns_read_part_force_hydro(casename{i}, time(tt));
                Hfx(tt,i) = mean(Fhx);
                Hfz(tt,i) = mean(Fhz);
                [Fix, Fiy, Fiz] = cgns_read_part_force_interaction(casename{i},time(tt));
                Ifx(tt,i) = mean(Fix);
                Ifz(tt,i) = mean(Fiz);
                [F1, F2, F3] = cgns_read_part_force_total(casename{i},time(tt));
                Fx(tt,i) = mean(F1);
                Fz(tt,i) = mean(F3);

                [U, V, W] = cgns_read_flow_vel(casename{i}, time(tt));
                [phase] = cgns_read_flow_phase(casename{i},time(tt));
                phase(phase~=-1) = 0;
                phase(phase==-1) = 1;
                phase = double(phase);
                [FUY, FUX, FUZ] = gradient(U,dy,dx,dz);
                [FVY, FVX, FVZ] = gradient(V,dy,dx,dz);
                [FWY, FWX, FWZ] = gradient(W,dy,dx,dz);

                %%dissipation of fluid is calculated as the velocity of particle is the solid body velocity
                epsilon_fluid(tt,i) = nu * sum(sum(sum(phase.*(2*FUX.*FUX  + 2*FVY.*FVY + 2*FWZ.*FWZ + FUY.*FUY + FVX.*FVX + 2*FUY.*FVX + ...
                                      FUZ.*FUZ + FWX.*FWX + 2*FUZ.*FWX + FVZ.*FVZ + FWY.*FWY + 2*FVZ.*FWY))))*dx*dy*dz;
                epsilon_part(tt,i) = -sum(Fhx.*u + Fhy.*v + Fhz.*w);

                %%potential eneregy of particle is calculated as (instantenous height - a),subtract the intial value
                PPE(tt,i) = sum((xe - x)*sin(theta_p{i}) + (z - zs)*cos(theta_p{i}))*(rho_p{i} - rho_f)*v_p*g - E_ref(i,1);
                PKE(tt,i) = 0.5 * sum(u.*u + v.*v + w.*w) * v_p * rho_p{i};
                Ex(tt,i) = 0.5 * sum(u.*u)*v_p * rho_p{i};
                Ez(tt,i) = 0.5 * sum(w.*w)*v_p * rho_p{i};

                U = U.*phase; V = V.*phase; W = W.*phase;
                FKE(tt,i) = 0.5 * sum(sum(sum((U.*U + V.*V + W.*W))))*m_f;
        end
end

Fhx = Hfx; Fhz = Hfz;
Fix = Ifx; Fiz = Ifz;

save('force_aspect.mat','time','norm_t','norm_f','Fhx','Fhz','Fix','Fiz','Fx','Fz');

for i = 1:4
        PPE(:,i) = PPE(:,i)/E0(i,1);
        PKE(:,i) = PKE(:,i)/E0(i,1);
        FKE(:,i) = FKE(:,i)/E0(i,1);
        Ex(:,i) = Ex(:,i)/E0(i,1);
        Ez(:,i) = Ez(:,i)/E0(i,1);
        epsilon_fluid(:,i) = epsilon_fluid(:,i)/E0(i,1);
        epsilon_part(:,i) = epsilon_part(:,i)/E0(i,1);
end

for tt = 1:Nt-1
        PPE_decay(tt,:) = (PPE(tt+1,:) - PPE(tt,:))/(time(1,tt+1) - time(1,tt));
        PKE_decay(tt,:) = (PKE(tt+1,:) - PKE(tt,:))/(time(1,tt+1) - time(1,tt));
        FKE_decay(tt,:) = (FKE(tt+1,:) - FKE(tt,:))/(time(1,tt+1) - time(1,tt));
end

save('energy_budget_aspect.mat','time','norm_t','PPE','PKE','Ex','Ez','FKE','PPE_decay','PKE_decay','FKE_decay','epsilon_fluid','epsilon_part');
