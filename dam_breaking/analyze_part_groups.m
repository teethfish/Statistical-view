clc;
close all;
clear all;

% Simulation name and time
casename{1} =
casename{2} = 
Ncase = length(casename);
bluePicpath = 
time = 
Nt = length(time);

% Groups in z direction: Pick up the highest Nm particles and store their sequence in index
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

% time loop: to find the mean velocity for these highest Nm particles
for tt = 1:Nt
    fprintf('time is %f\n',time(tt))
    for j = 1:Ncase
        [u, v, w] = cgns_read_part_vel(casename{j},time(tt));
        avg_Nm_vel_x(tt,j) = u.*index/Nm;
        avg_Nm_vel_z(tt,j) = w.*index/Nm;
    end
end    
        

% Groups in x directions
[x, y, z] = cgns_read_part_position(casename{1},0);
[~, index] = sort(x);
right_index = index(1:80);
middle_index = index(81:160);
left_index = index(161:240);

right_ke_u = zeros(Nt,1);
middle_ke_u = zeros(Nt,1);
left_ke_u = zeros(Nt,1);
right_ke_w = zeros(Nt,1);
middle_ke_w = zeros(Nt,1);
left_ke_w = zeros(Nt,1);

% time loop: to find the mean kinetic energy of three different groups
for tt = 1:Nt
    fprintf('time is %f\n',time(tt));
    for j = 1:Ncase
        [u, v, w] = cgns_read_part_vel(casename{j},time(tt));
        for i = 1:80
            right_ke_u(tt,j) = right_ke_u(tt,1) + u(right_index(i))*u(right_index(i));
            middle_ke_u(tt,j) = middle_ke_u(tt,1) + u(middle_index(i))*u(middle_index(i));
            left_ke_u(tt,j) = left_ke_u(tt,1) + u(left_index(i))*u(left_index(i));

            right_ke_w(tt,j) = right_ke_w(tt,1) + w(right_index(i))*w(right_index(i));
            middle_ke_w(tt,j) = middle_ke_w(tt,1) + w(middle_index(i))*w(middle_index(i));
            left_ke_w(tt,j) = left_ke_w(tt,1) + w(left_index(i))*w(left_index(i));
        end
    end
end


