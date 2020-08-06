clc
clear all

% diffusivity 
k = 0.5;
% length of x
L_x = 10;
% simlation time
t = 12;
% the number of grids
grid_t = 200;
grid_x = 20;

delta_x = L_x / (grid_x + 1);
delta_t = t / (grid_t + 1);

A = diag(-4*ones(1,grid_x^2)) + diag(ones(1,grid_x^2-1),1) + diag(ones(1,grid_x^2-1),-1);
A = A + diag(ones(1,grid_x^2-20),20) + diag(ones(1,grid_x^2-20),-20);
A = k*A*(1/delta_x)^2;

% Dirichlet boundry conditions
b = zeros(grid_x^2,1);
b = (1/delta_x)^2*b*k;

% define initial state
x = 0:delta_x:L_x; 
y = x;
[X,Y] = meshgrid(x,y);
R = sqrt((X- L_x/2).^2 + (Y-L_x/2).^2);
Z = 100*sin(R)./R;

% solving
data = Z;
for i = 1:grid_t
    new_data = zeros(grid_x+2,grid_x+2);
    data(:,:,i+1) = new_data;
    prev_v = reshape(data(2:grid_x+1,2:grid_x+1,i)',[grid_x^2,1]);
    v = linsolve((eye(grid_x^2)/delta_t-A),prev_v/delta_t+b)';
    data(2:grid_x+1,2:grid_x+1,i+1) = reshape(v,[grid_x,grid_x])';
end

% visualization: animation
taxis = linspace(0,t,grid_t+1);
for m = 1:(grid_t+1)
    s = mesh(X,Y,data(:,:,m),'FaceAlpha','0.5');
    s.FaceColor = 'flat';
    title(sprintf('Time elapsed: %0.2f sec', taxis(m)));
    xlim([0 L_x]); 
    ylim([0 L_x]);
    zlim([-100,100]);
    drawnow;   
end