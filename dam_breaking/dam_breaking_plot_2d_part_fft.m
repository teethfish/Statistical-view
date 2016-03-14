clear all; close all; clc;

casename = '/home-4/ywang221@jhu.edu/scratch/dam_breaking/sedimentation_N240/as1.25';
bluePicPath = '/home-4/ywang221@jhu.edu/bluePic/dam_breaking/N240_material_property';
time = 1:1:26;
Nt = length(time);

% particle radius and half length of domain in three directions
a = 1; 
Hx = 30; Hy = 4.5; Hz = 15;

% expansion order and reconstructed intervals
M = 30; NUM_x = 60;
L = 15; NUM_z = 30;
xx = linspace(-Hx,Hx,NUM_x);xx = xx';
zz = linspace(-Hz,Hz,NUM_z);zz = zz';

for tt = 1 : Nt
        fprintf('time is %f\n',time(tt));

        [u, v, w] = cgns_read_part_vel(casename, time(tt));
        [x, y, z] = cgns_read_part_position(casename, time(tt));

        % kinetic energy in two directions
        Ex = 0.5 * u.* u;
        Ez = 0.5 * w.* w;

        Np = length(u);
        density = ones(Np,1);

        % calculate the fourier coefficents
        [nq(:,:,tt)] = dam_breaking_calculate_part_2d_fft(Ez, x, z, Hx, Hz, M, L);
        [nq_density(:,:,tt)] = dam_breaking_calculate_part_2d_fft(density, x, z, Hx, Hz, M, L);

        % reconstruct the kinetic energy distribution
        [NQ(:,:,tt)] = dam_breaking_reconstruct_fft_2d(xx, zz, M, L, nq(:,:,tt), Hx, Hz, a, 0);
        [NQ_density(:,:,tt), beta(:,:,tt)] = dam_breaking_reconstruct_fft_2d(xx, zz, M, L, nq_density(:,:,tt), Hx, Hz, a, 0);
        %beta(:,:,tt) = beta(:,:,tt)/2/Hy;
        %NQ(:,:,tt) = NQ(:,:,tt)./NQ_density(:,:,tt);

        handle = figure;set(gca,'fontsize',16);
        contourf(zz+1,xx,NQ(:,:,tt));
        caxis([0, 0.5]);colorbar;
        str = sprintf('vertical kinetic energy, t = %.2f',time(tt));title(str);
        view([90 -90]);xlabel('z');ylabel('x')
        filename=[bluePicPath,'/','vertical_kenetic_energy_M=',num2str(M),'L=',num2str(L),'t=',num2str(time(tt)),'.eps'];
        saveas(handle,filename,'epsc');
end
