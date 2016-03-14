%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script is used to get the two-particle distribution function
%
%INPUTS:
%   Nfile:  the total number of realization, each one stored in a file with
%   Hx:     length of domain in x direction
%   Hy:     length of domain in y direction
%   Hz:     length of domain in z direction
%   d:      diameter of particle
%   rgmax:  we're intersted in the distance within rgmax*diameter
%   mg:     divide the distance into "mg" parts
%Parameters:
%   NP:     the total number of particles
%   vf:     volume fraction
%   ddcst:  normalization of g(r)
%   ddr:    normalization of g(r)
%RETURNS:
%   r:      the normalized(/d) distance in radial direction
%   g:      the radius distribution function we calculated
%Uses: function [plot_r, p_y] = percus_yevik(v_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all

Nfile = 1;
for i=1:Nfile
    filename{i} = ['/Users/yayun/Documents/MATLAB/dataprocessing/seeder/test_num_of_real/N=1/myrandom_' num2str(i-1,'%d') '.txt'];
    %filename{i} = ['/Users/yayun/Documents/MATLAB/dataprocessing/seeder/NP=27/myrandom_' num2str(i+1000,'%d') '.txt'];
end

A = importdata(filename{1});
x = A(:,1);
NP = length(x)

d = 1; %diameter of particles
Hx = 5;
Hy = 5;
Hz = 5;
V = Hx*Hy*Hz;

%get the solution of percus-Yevik
v_f = NP*pi*d*d*d/(6*Hx*Hy*Hz);
[plot_r, p_y] = percus_yevik(v_f);

H = min([Hx, Hy, Hz]);
rgmax = 5;
srgmax = rgmax*d;
if srgmax >= H/2
  srgmax = H/2;
end
rgmax=srgmax/d;

%get normalization of g(r)
mg=20;
ddcst=4.0*pi*(NP-1.0)*NP*d*d*d/3.0/V;
dgr=(rgmax-1.0)/mg;
r=zeros(1,mg);
ddr=zeros(1,mg);
f=zeros(1,mg);
g=zeros(1,mg);
for n=1:mg
  rl=1.0+n*dgr;
  rs=rl-dgr;
  r(n)=sqrt((rl*rl+rs*rs)/2.0);
  ddr(n)=rl*rl*rl-rs*rs*rs;
end
erow=ones(1,NP);
ecol=ones(NP,1);


for n=1:Nfile
    A = importdata(filename{n});
    x = A(:,1);
    y = A(:,2);
    z = A(:,3);
    rx=x*erow-ecol*x';
    ry=y*erow-ecol*y';
    rz=z*erow-ecol*z';
    rx=rx-Hx*fix(2.0*rx/Hx);
    ry=ry-Hy*fix(2.0*ry/Hy);
    rz=rz-Hz*fix(2.0*rz/Hz);
    rr2=rx.^2+ry.^2+rz.^2;
    ix=fix((sqrt(rr2/d/d)-1.0)/dgr)+1;
    for i=1:NP
        for j=1:NP
            if ix(i,j) >= 1 && ix(i,j) <= mg
                f(ix(i,j))=f(ix(i,j))+1.0;
            end
        end
    end
end
for jj=1:mg
  g(jj)=f(jj)/ddcst/ddr(jj)/Nfile;
end

g
figure
set(gca,'fontsize',16)
plot(r,g,'k*','linewidth',2)
hold on
%plot(interval(1:M),rdf(1),'k','linewidth',2)
plot(plot_r, p_y,'r','linewidth',2)
legend('simulation','Percus-Yevik')
title(['Num of realization=' num2str(Nfile)]);






