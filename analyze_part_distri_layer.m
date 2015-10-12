%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script will generally separate the particle into groups due to their inital position
%and track the mean velocity and position in each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all;

casename = 
time = 0:0.1:30;
Nt = length(time);
savename = '_as=1';
savedata = 'part_layer.mat';

[xs, ~, ~, ~, ~, ~, ~, ~, zs] = calculate_domain_info(casename);

%% generate particle group, particles are saves from bottom to top
%% part_number(M,Ncase) will record how many particles in each layer, part(tmp, j, i) will store particle's index in each layer
%% M is the particle layers
M = 8; 
part_number = zeros(M,1);
[x, y, z] = cgns_read_part_position(casename, 0.0);
Np = length(z);
z_min = min(z); z_max = max(z) + 0.01;
delta_z = (z_max - z_min)/M;
for j = 1:M
        for p = 1:Np
                if z(p) >= z_min + (j-1)*delta_z && z(p)< z_min + j*delta_z
                        part_number(j,1) = part_number(j,1) + 1;
                        tmp = part_number(j,1);
                        part(tmp,j) = p;
                end
        end
end


%accumulate statistics
w_mean = zeros(Nt, M);
z_mean = zeros(Nt, M);

for tt = 1:Nt
        fprintf('time is %f\n',time(tt));
                
        [u, v, w] = cgns_read_part_vel(casename, time(tt));
        [x, y, z] = cgns_read_part_position(casename, time(tt));
        for j = 1:M
                for mm = 1:part_number(j,1)
                        tmp = part(mm,j);
                        %scatter_w(mm,j,tt,i) = w(tmp); %this is for scatter plot, snapshots of w versus x position
                        %scatter_x(mm,j,tt,i) = x(tmp);
                        tmp = part(mm,j);
                        u_mean(tt,j) = u_mean(tt,j) + u(tmp); %x velocity in each layer 
                        w_mean(tt,j) = w_mean(tt,j) + w(tmp); %z velocity in each layer
                        x_mean(tt,j) = x_mean(tt,j) + x(tmp); %wave front position in each layer
                        z_mean(tt,j) = z_mean(tt,j) + z(tmp); %decrease 
                end
                u_mean(tt,j) = u_mean(tt,j)/part_number(j,1);
                w_mean(tt,j) = w_mean(tt,j)/part_number(j,1);
                x_mean(tt,j) = x_mean(tt,j)/part_number(j,1);
                z_mean(tt,j) = z_mean(tt,j)/part_number(j,1);
        end
end

%rescale x position and z position by subtract the reference bottom to left wall
z_mean = z_mean - zs;
x_mean = x_mean - xs;

eval(sprintf('time%s = time;',savename));
eval(sprintf('u_mean%s = u_mean;',savename));
eval(sprintf('w_mean%s = w_mean;',savename));
eval(sprintf('x_mean%s = x_mean;',savename));
eval(sprintf('z_mean%s = z_mean;',savename));

save(savedata,'time_*','u_mean_*','w_mean_*','x_mean_*','z_mean_*');

