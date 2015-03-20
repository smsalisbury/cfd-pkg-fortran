close all;

% Run fortran
! a.out

% Get variables
u = dlmread('outputs/u.out');
v = dlmread('outputs/v.out');
velocity = (u.^2+v.^2).^0.5;
cont = dlmread('outputs/continuity.out');
axes = dlmread('outputs/P_axes.out');
x = axes(:,1);
y = axes(:,2);

% Make plots
figure;
semilogy(cont(:,1),cont(:,2));
figure;
quiver(x,y,u,v,2);
hold on;
contour(x,y,velocity);
