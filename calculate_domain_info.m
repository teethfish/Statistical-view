function [Xs, Xe, Lx, dx, Ys, Ye, Ly, dy, Zs, Ze, Lz, dz] = calculate_domain_info(casename)

[X, Y, Z] = cgns_read_grid(casename);
[Xn, Yn, Zn] = size(X);
Xn = Xn - 1;
Yn = Yn - 1;
Zn = Zn - 1;

Xs = X(1,1,1);
Xe = X(end,1,1);
Lx = Xe - Xs;
dx = Lx/Xn;

Ys = Y(1,1,1);
Ye = Y(1,end,1);
Ly = Ye - Ys;
dy = Ly/Yn;

Zs = Z(1,1,1);
Ze = Z(1,1,end);
Lz = Ze - Zs;
dz = Lz/Zn;

end
