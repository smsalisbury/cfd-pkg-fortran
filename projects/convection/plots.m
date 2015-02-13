close all;

% Get variables
phi = dlmread('output/phi.out');
axes = dlmread('output/axes.out');
conv = dlmread('output/conv.out');
v = [-10 10 20 30 40 50 60 70 80 90];
cross = [];
for i=1:size(phi,1)
    cross(i) = phi(i,size(phi,1)-i+1);
end
cross = transpose(cross)

% Make plots
figure;
contourf(axes(:,1),axes(:,2),phi,v);
figure;
semilogy(conv(:,1),conv(:,2));
figure;
plot(axes(:,1),cross);