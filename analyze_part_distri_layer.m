%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script will generally separate the particle into groups due to their inital position
%and tract the mean velocity and position in each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all;

casename{1} = 
casename{2} = 
Ncase = length(casename);
[xs, ~, ~, ~, ~, ~, ~, ~, zs] = calculate_domain_info(casename{1});

time = 0:0.1:30;
Nt = length(time);

%% generate particle group, particles are saves from bottom to top
%% part_number(M,Ncase) will record how many particles in each layer, part(tmp, j, i) will store particle's index in each layer
%% M is the particle layers
M = 8; 
part_number = zeros(M,Ncase);
for i = 1:Ncase
        [x, y, z] = cgns_read_part_position(casename{i}, 0.0);
        Np = length(z);
        z_min = min(z); z_max = max(z) + 0.01;
        delta_z = (z_max - z_min)/M;
        for j = 1:M
                for p = 1:Np
                        if z(p) >= z_min + (j-1)*delta_z && z(p)< z_min + j*delta_z
                                part_number(j,i) = part_number(j,i) + 1;
                                tmp = part_number(j,i);
                                part(tmp,j,i) = p;
                        end
                end
        end
end

%accumulate statistics
w_mean = zeros(Nt, M, Ncase);
z_mean = zeros(Nt, M, Ncase);

for tt = 1:Nt
        fprintf('time is %f\n',time(tt));
        for i = 1:1
                [u, v, w] = cgns_read_part_vel(casename{i}, time(tt));
                [x, y, z] = cgns_read_part_position(casename{i}, time(tt));
                for j = 1:M
                        for mm = 1:part_number(j,i)
                                tmp = part(mm,j,i);
                                %scatter_w(mm,j,tt,i) = w(tmp); %this is for scatter plot, snapshots of w versus x position
                                %scatter_x(mm,j,tt,i) = x(tmp);
                                tmp = part(mm,j,i);
                                u_mean(tt,j,i) = u_mean(tt,j,i) + u(tmp); %x velocity in each layer 
                                w_mean(tt,j,i) = w_mean(tt,j,i) + w(tmp); %z velocity in each layer
                                x_mean(tt,j,i) = x_mean(tt,j,i) + x(tmp); %wave front position in each layer
                                z_mean(tt,j,i) = z_mean(tt,j,i) + z(tmp); %decrease 
                        end
                        u_mean(tt,j,i) = u_mean(tt,j,i)/part_number(j,i);
                        w_mean(tt,j,i) = w_mean(tt,j,i)/part_number(j,i);
                        x_mean(tt,j,i) = x_mean(tt,j,i)/part_number(j,i);
                        z_mean(tt,j,i) = z_mean(tt,j,i)/part_number(j,i);
                end
        end
end

%rescale x position and z position by subtract the reference bottom to left wall
z_mean = z_mean - zs;
x_mean = x_mean - xs;

save('part_layer.mat','time','u_mean','w_mean','x_mean','z_mean');

