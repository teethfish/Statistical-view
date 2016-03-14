%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This main script is used to calculate g_l0 of radial distribution function
%Notice:g_00 is the g(r), averaged over the entire theta
%Input:
%       d:          diameter
%       Hx,Hy,Hz:   domain length in three direction
%       xs,ys,zs:   the start point in three direction
%       mg:         divide the distance we're intereted into small parts
%Output:
%       r:          for plot purpose
%       g:          radial distribution function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; close all; clear all

bluePicPath = '/Users/yayun/Documents/MATLAB/dataprocessing/seeder/bias_move'; %change the path
Nfile = 19;  %number of realization
for i=1:Nfile
    filename{i} = ['/Users/yayun/Documents/MATLAB/dataprocessing/seeder/bias_move/myrandom_' num2str(i-1,'%d') '.txt'];
end

A = importdata(filename{1});
x = A(:,1);
NP = length(x);

d = 1; %diameter of particles
Hx = 5;Hy = 5;Hz = 5; 
xs = 0;ys = 0;zs = 0;
V = Hx*Hy*Hz;
v_f = NP*pi*d*d*d/(6*Hx*Hy*Hz);
[plot_r, p_y] = percus_yevik(v_f);%get the solution of percus-Yevik

%Decide the maximum domain can be considered(periodicity consideration)
H = min([Hx, Hy, Hz]);
rgmax = 5;
srgmax = rgmax*d;
if srgmax >= H/2
  srgmax = H/2;
end
rgmax=srgmax/d;

%calculate r(n) for drawing
mg = 40;                  %divide the region into mg parts
dgr = (rgmax - 1.0)/mg;
r = zeros(1,mg);
for n=1:mg
  rl = 1.0 + n*dgr;
  rs = rl - dgr;
  r(n) = sqrt((rl*rl + rs*rs)/2.0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate g_l0(r),only a function of r
gl = zeros(mg,9); %in three direction, each one has three components x,y,z,we have gl_0, gl_2, gl_4
for nn=1:Nfile                  %average the results for N realization   
    A = importdata(filename{nn});
    x = A(:,1);
    y = A(:,2);
    z = A(:,3);
    for ll = 1:3        %Here ll can be changed to get different order of g(r,theta)
        l = 2*(ll-1);
        index = (ll-1)*3;
        [g_x,g_y,g_z] = function_g_l0(NP, l, d, dgr, x, y, z, mg, Hx, Hy, Hz, xs, ys, zs);
        gl(:,index+1) = gl(:,index+1) + g_x*sqrt((2*l+1)/(4*pi));
        gl(:,index+2) = gl(:,index+2) + g_y*sqrt((2*l+1)/(4*pi));
        gl(:,index+3) = gl(:,index+3) + g_z*sqrt((2*l+1)/(4*pi));
    end
end
gl = gl/Nfile;

%plot gl_0,gl_2 and gl_4
handle=figure;
set(gca,'fontsize',16)
plot(r,gl(:,1),'r*','linewidth',2)
hold on
plot(plot_r, p_y,'k','linewidth',2)
legend('Polar-z','Percus-Yevik')
title('g-00 realization=1');    %change the title
filename=[bluePicPath,'/','g00_polar_z(1)','.png'];
saveas(handle,filename);

handle=figure;
set(gca,'fontsize',16)
plot(r,gl(:,4),'k-','linewidth',2)
legend('Polar-z')
title('g-20 realization=1');
filename=[bluePicPath,'/','g20_polar_z(1)','.png'];
saveas(handle,filename);

handle=figure;
set(gca,'fontsize',16)
plot(r,gl(:,7),'k-','linewidth',2)
legend('Polar-z')
title('g-40 realization=1');
filename=[bluePicPath,'/','g40_polar_z(1)','.png'];
saveas(handle,filename);
% for ll=1:3
%     index = (ll-1)*3;  
%     figure
%     set(gca,'fontsize',16)
%     plot(r,gl(:,1+index),'r','linewidth',2)
%     legend('Polar-x')
%     title('realization=1');
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %calculate g(r,theta), also consider the theta
% M = 30;
% theta = linspace(0,pi/2,M);
% g = zeros(mg,1);
% sum_gr = zeros(mg,M);
% for nn=1:Nfile
%     A = importdata(filename{nn});
%     x = A(:,1);
%     y = A(:,2);
%     z = A(:,3);
%     [tmp_g] = function_g_l0(NP, 0, d, dgr, x, y, z, mg, Hx, Hy, Hz, xs, ys, zs);
%     g = g + tmp_g*sqrt(1/(4*pi));
%     for ll = 1:4               %Here ll can be changed to get different order of g(r,theta)
%         l = 2*(ll-1);
%         [g_x,g_y,g_z] = function_g_l0(NP, l, d, dgr, x, y, z, mg, Hx, Hy, Hz, xs, ys, zs);
%         for mm=1:M
%             mu = cos(theta(mm)); 
%             sum_gr(:,mm) = sum_gr(:,mm) + sqrt((2*l+1)/(4*pi))*legendre(l,mu)*g_y;  %check which axis is the polar axis
%         end
%     end
% end
% g = g/Nfile;    
% sum_gr = sum_gr/Nfile;
% 
% %divide sum_gr by g(r) to get only the theta dependence
% for i=1:mg
%     sum_gr(i,:) = sum_gr(i,:)/g(i);
% end
% 
% %CONTOUR plot of g(r,theta) function
% th = linspace(pi/2,0,M);
% [TH,R] = meshgrid(th,r);
% [X,Y] = pol2cart(TH,R);
% handle=figure;
% set(gca,'FontSize',16);
% surf(X,Y,sum_gr)
% title('Polar-y')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








