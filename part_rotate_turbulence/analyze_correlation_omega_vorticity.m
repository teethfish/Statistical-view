clc; close all; clear all;

%Simulation name and time(time should be starts from the stationary state)
casename = '/home/yawang/bluebottle_0.1/cases/turb_rotation_convergence_test/dowm_3a_up5a';
t = 0:0.005:2.5;
dt = t(1,2) - t(1,1);
Nt = length(t);

%Flow property: turbulence forcing and inlet velocity; particle radius
nu = 1.0; turbA = 7;
a = 1;
W_mean = 75;

%Choose planes for considering w(x=x_p) = x(x=0,t-Ld/W_mean)
Ld = [3, 2];
t_lag = Ld*a/W_mean;

%Domain information
[~, ~, Lx, dx, Nx, ~, ~, Ly, dy, Ny, ~, ~, Lz, dz, Nz] = calculate_domain_info_turb(casename);
zn = 105 - round(Ld*a/dz);
Mz = length(zn);

%Generate statistics
for tt = 1:Nt
	fprintf('time is %f\n',t(tt));	

	%Check the time for calculating vorticity
	time = t(tt) - t_lag;

	%find the data points for vorticity
	array = abs(time - t);
	[~, Mt] = min(array);

	%Calculating vorticity on time (t-t_lag)
	[u, v, w] = cgns_read_flow_vel(casename, t(Mt));

	% substract the mean vertical velocity
	w = w - W_mean; 

	% Calculate fluid vorticity on plane
	[wx(tt,:,:,:), wy(tt,:,:,:), wz(tt,:,:,:)] = calculate_fluid_vorticity_on_plane(u, v, w, dx, dy, dz, zn);

	%Get particle scalar
	[Ox(tt), Oy(tt), Oz(tt)] = cgns_read_part_omega(casename, t(tt));
	[Couple_x(tt), Couple_y(tt), Couple_z(tt)] = cgns_read_part_moment(casename, t(tt));
end

%Do the correlation for each nodes on the plane:for each point do a correlation of (OmegaVorticiey)
for zz = 1:Mz
	for i = 1:Nx
		for j = 1:Ny
			[S, Cow_x(:,i,j)] = calculate_1d_correlation(Ox, wx(:,i,j,zz));
			[~, Cow_y(:,i,j)] = calculate_1d_correlation(Oy, wy(:,i,j,zz));
			[~, Cow_z(:,i,j)] = calculate_1d_correlation(Oz, wz(:,i,j,zz));

			[~, Clw_x(:,i,j)] = calculate_1d_correlation(Couple_x, wx(:,i,j,zz));
			[~, Clw_y(:,i,j)] = calculate_1d_correlation(Couple_y, wy(:,i,j,zz));
			[~, Clw_z(:,i,j)] = calculate_1d_correlation(Couple_z, wz(:,i,j,zz));	
		end
	end
end


Cow_x_avg = mean(mean(Cow_x,2),2);
Cow_y_avg = mean(mean(Cow_y,2),2);
Cow_z_avg = mean(mean(Cow_z,2),2);

Clw_x_avg = mean(mean(Clw_x,2),2);
Clw_y_avg = mean(mean(Clw_y,2),2);
Clw_z_avg = mean(mean(Clw_z,2),2);
