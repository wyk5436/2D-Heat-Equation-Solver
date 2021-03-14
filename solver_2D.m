clc
clear all

% diffusivity 
k = 0.5;
% length of x
L_x = 10;
% simulation time
t = 6;
% the number of grids
grid_t = 200;
grid_x = 100;

delta_x = L_x / (grid_x + 1);
delta_t = t / (grid_t + 1);

e = ones(grid_x^2,1);
A = spdiags([e -4*e e],-1:1,grid_x^2,grid_x^2);
A = A + spdiags(e,grid_x,grid_x^2,grid_x^2) + spdiags(e,-grid_x,grid_x^2,grid_x^2);
I = spdiags(e,0,grid_x^2,grid_x^2);
for i = 1:grid_x-1
    A(grid_x*i,grid_x*i+1) = 0;
    A(grid_x*i+1,grid_x*i) = 0;
end
A = A*((k*delta_t)/(2*delta_x^2));

% Dirichlet boundry conditions
b = zeros(grid_x^2,1);
b = (delta_t/(2*delta_x^2))*b*2*k;

% define initial state
x = 0:delta_x:L_x; 
y = x;
[X,Y] = meshgrid(x,y);
R = sqrt((X- L_x/2).^2 + (Y-L_x/2).^2);
Z = 100*sin(R)./R;
% enforcing boundary condition
Z(1,:) = 0;
Z(grid_x+2,:) = 0;
Z(:,1) = 0;
Z(:,grid_x+2) = 0;

% solving
data = Z;
for i = 1:grid_t
    new_data = zeros(grid_x+2,grid_x+2);
    data(:,:,i+1) = new_data;
    prev_v = reshape(data(2:grid_x+1,2:grid_x+1,i)',[grid_x^2,1]);
    RHS = (A+I)*prev_v + b;
    v = (I-A)\RHS;
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
