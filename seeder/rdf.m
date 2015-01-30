%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script is used to calculate the radial distribution function and
%compare it with the Percus-Yevik solution.
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

Nfile = 1;
for i=1:Nfile
    filename{i} = ['/Users/yayun/Documents/MATLAB/dataprocessing/seeder/test_num_of_real/N=1/myrandom_' num2str(i-1,'%d') '.txt'];
end

A = importdata(filename{1});
x = A(:,1);
NP = length(x);

d = 1; %diameter of particles
Hx = 5;Hy = 5;Hz = 5; 
xs = 0;ys = 0;zs = 0;
V = Hx*Hy*Hz;

%get the solution of percus-Yevik
v_f = NP*pi*d*d*d/(6*Hx*Hy*Hz);
[plot_r, p_y] = percus_yevik(v_f);

%Decide the maximum domain can be considered(periodicity consideration)
H = min([Hx, Hy, Hz]);
rgmax = 5;
srgmax = rgmax*d;
if srgmax >= H/2
  srgmax = H/2;
end
rgmax=srgmax/d;

%calculate r(n) for drawing
mg=20;
dgr=(rgmax-1.0)/mg;
r=zeros(1,mg);
for n=1:mg
  rl=1.0+n*dgr;
  rs=rl-dgr;
  r(n)=sqrt((rl*rl+rs*rs)/2.0);
end

%calculate g(r), consider the periodicity
g = zeros(mg,1);
for nn=1:Nfile
    A = importdata(filename{nn});
    x = A(:,1);
    y = A(:,2);
    z = A(:,3);
    for i=1:NP
        for j=1:NP
            dx=x(j)-x(i);dz=z(j)-z(i);dy=y(j)-y(i);
            if abs(x(j)-x(i))>Hx/2
                dx=Hx-abs(x(j)- x(i));
            end
            if abs(y(j)-y(i))>Hy/2
                dy=Hy-abs(y(j)- y(i));
            end
            if abs(z(j)-z(i))>Hz/2
                dz=Hz-abs(z(j)- z(i));
            end 
            d_two = dx*dx+dy*dy+dz*dz;
            tmp_d = fix((sqrt(d_two)/d - 1)/dgr)+1;
            if tmp_d > 0 && tmp_d <= mg
              g(tmp_d) = g(tmp_d) + 1/d_two;  
            end
        end
    end
            
end
%Normalize g(r)
g = g*V/(4*pi*NP*(NP-1)*dgr)/Nfile;

g
figure
set(gca,'fontsize',16)
plot(r,g,'k*','linewidth',2)
hold on
%plot(interval(1:M),rdf(1),'k','linewidth',2)
plot(plot_r, p_y,'r','linewidth',2)
legend('simulation','Percus-Yevik')
title(['Num of realization=' num2str(Nfile)]);








