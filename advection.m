function [Pav] = advection(A,B,L,Ncells,t_final)
% This function solves 1D transport equation using finite volume method. 
%   Input arguments
%       A,B     Constants in Gaussian function (as initial condition)
%       L       The length of interest
%       Ncells  Number of cells
%       dt      Time step
%       t_final Final time
%
%   Output arguments
%       P    Array of average densities

% Avearge velocity
u = 1;

% Set the cells
dx = L/Ncells;
Pav = zeros(1,Ncells);
dPav = zeros(1,Ncells);

% CFL condition
dt = u*dx;
x = 0:dx:L;
Ncells = L/dx;
xav = 0.5*(x(1:Ncells)+x(2:Ncells+1));

% Calculate number of steps
Nsteps = round(t_final/dt); 

% Initial condition
for j = 1:Ncells

    Pav(j) = integral(@(x) exp(-A*(x-B).^2),x(j),x(j+1))/dx;


end

P1=Pav;
t = zeros(Nsteps,1);
% Calculate the constants
C1 = u/dx;

writerObj = VideoWriter('Advection');
writerObj.FrameRate = 20;

writerObj.Quality = 100;
open(writerObj);

% Main loop
for i = 1:Nsteps
    
    clf
    for j = 1:Ncells
    
    % Forward method
    % dPav(j) = -C1*(Pav(j+1)-Pav(j));
    
    % Backward method
    dPav(j) = -C1*(Pav(j)-Pav(mod(j-2,Ncells)+1));
   
    % Central difference method
    % dPav(j) = -2*C1*(Pav(j+1)-Pav(j-1));
     
    end
 
    
    Pav = Pav + dPav*dt;
    
    t(i+1) = t(i) +dt;
    
    plot(xav,Pav);
    
    xlim([0 L]);
   
    frame = getframe(gcf);
    writeVideo(writerObj,frame);


end
close(writerObj);
end