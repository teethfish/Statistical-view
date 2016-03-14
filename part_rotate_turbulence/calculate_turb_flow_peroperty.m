clc; close all; clear all;

%Simulation name and time
casename = '/home/yawang/bluebottle_0.1/cases/turb_rotation_convergence_test/dowm_3a_up5a';
t = 0.23:0.01:0.46;
dt = t(1,2) - t(1,1);
Nt = length(t);

%Flow property: turbulence forcing and inlet velocity
nu = 1.0; turbA = 7;
W_mean = 75;

%Choose planes for results
zn = [1,8,15,23,30,38,45,53,59,91,98,105,113];
Mz = length(zn);

%Domain information
[~, ~, Lx, dx, Nx, ~, ~, Ly, dy, Ny, ~, ~, Lz, dz, Nz] = calculate_domain_info_turb(casename);

%Generate statistics
for tt = 1:Nt
	fprintf('time is %f\n',t(tt));
	[u, v, w] = cgns_read_flow_vel(casename, t(tt));

	% substract the mean vertical velocity
	w = w - W_mean; 

	% calculate dissipation rate on each plane and overall dissipation
	[epsilon(tt,:), epsilon_overall(tt,1)] = calculate_fluid_dissipation_on_plane(u, v, w, dx, dy, dz, nu, zn);
	
	
	for zz = 1:Mz
		fprintf('calculating plane %d\n',zz);
		U = u(:,:,zn(zz)); V = v(:,:,zn(zz)); W = w(:,:,zn(zz));
		u_mean(tt,zz) = mean(mean(U));v_mean(tt,zz) = mean(mean(V));w_mean(tt,zz) = mean(mean(W));
        u_fluc = U - u_mean(tt,zz); v_fluc = V - v_mean(tt,zz); w_fluc = W - w_mean(tt,zz);
        u_rms(tt,zz) = sqrt(mean(mean(u_fluc.*u_fluc)));
        v_rms(tt,zz) = sqrt(mean(mean(v_fluc.*v_fluc)));
        w_rms(tt,zz) = sqrt(mean(mean(w_fluc.*w_fluc)));
        tkn(tt,zz) = 0.5 * mean(mean((u_fluc.*u_fluc + v_fluc.*v_fluc + w_fluc.*w_fluc)));
        kolm_l(tt,zz) = (nu^3/epsilon(tt,zz))^(0.25);
    end
end

turnover_t = tkn./epsilon;
taylor_l = sqrt(15*nu*(w_rms.*w_rms)./epsilon);
int_l = 0.9*(w_rms.*w_rms).*w_rms./epsilon;
Re_g = w_rms.*taylor_l/nu;
