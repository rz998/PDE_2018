function [T] = heat(BC,IC,L,Nnodes,t_final)
% This function solves 1D heat equation using finite difference method. 
%   Input arguments
%       BC      Boundary condition
%       IC      Initial condition
%       L       The length of interest
%       Nnodes  Number of nodes
%       dt      Time step
%       t_final Final time
%
%   Output arguments
%       T    Array of temperatures 

alpha = 0.0005;

% Set the cells
dx = L/(Nnodes-1);
x = 0:dx:L;
dt = 0.5*dx^2/2/alpha;
Nsteps = round(t_final/dt);
t = zeros(Nsteps,1);

% Initial condition
T = ones(Nnodes-2,1)*IC;
T = [BC(1);T;BC(2)];

dT = zeros(Nnodes,1);

% Calculate the constants
C1 = -2*alpha/dx^2;
C2 = alpha/dx^2;

% Analytical solution
T_a = x*(BC(2)-BC(1))/L + BC(1);


writerObj = VideoWriter('Heat');
writerObj.FrameRate = 20;

writerObj.Quality = 100;
open(writerObj);

% Main loop
for i = 1:Nsteps
    
    clf 
    
    for j = 2:Nnodes-1
    
    dT(j) = C2*(T(j-1)+T(j+1))+C1*T(j);
    
    end
    
    % Evolve the temperature 
    T = T + dT*dt;
    
    % Advance in time
    t(i+1) = t(i) + dt;
    
    plot(x,T_a);
    hold on
    plot(x,T);
    
    xlim([0 L]);
   
    frame = getframe(gcf);
    writeVideo(writerObj,frame);

    
end
close(writerObj);
end
