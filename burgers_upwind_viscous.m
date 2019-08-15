function [U] = burgers_upwind_viscous(L,Ncells,Nsteps)
% This function solves 1D Burgers' equation using the upwind scheme
% under viscous condtions.
%   Input arguments
%       L       The length of interest
%       Ncells  Number of cells
%       Nsteps  Number of time steps
%
%   Output arguments
%       U       Velocity profile

% Viscosity
v = 0.00002;

% Set the cells
dx = L/Ncells;
x = 0:dx:L;
xav = 0.5*(x(1:Ncells)+x(2:Ncells+1));
U = zeros(1,Ncells);
t = zeros(5000,1);


% Calculate number of steps
Nsteps = 5000; 

% Initial condition
for j = 1:Ncells

    U(j) = integral(@(x) exp(-5*(x-1).^2),x(j),x(j+1))/dx;

end

% Video writer
writerObj = VideoWriter('Burgers_upwind_viscous');
writerObj.FrameRate = 20;

writerObj.Quality = 100;
open(writerObj);

% Main loop
for i = 1:Nsteps
    
    % Find the gradient 
    
    [dU2] = upwind(U);
    
    % Compute the viscous term
    for j = 1:Ncells
        
        d2U = U(mod(j,Ncells)+1) - 2*U(j) + U(mod(j-2,Ncells)+1);
    
    end
    
    % CFL condition (Cmax = 1)

    dt = dx/max(abs(U));

    % Evolve with respect to time

    U = U + dt*((0.5*-dU2/dx) + v*d2U/dx^2);
    
    plot(xav,U);

    xlim([0 L]);
 
    frame = getframe(gcf);
    writeVideo(writerObj,frame)



end

    function [dU2] = upwind(U)
    % This function uses the minmod limiter to find the gradient 
    %   Input arguments
    %       
    %
    %   Output arguments
    %       dU2          Selected gradient 
    
    dU2 = zeros(1,Ncells);
    
    for k = 1:Ncells
        
        if U(k) > 0
            
            dU2(k) = U(k)^2-U(mod((k-2),Ncells)+1)^2;
        
        elseif U(k) < 0
            
            dU2(k) = U(mod((k-2),Ncells)+1)^2-U(k)^2;
            
        else
            
            dU2(k) = 0;
        
        end
     
    end
    
    end    
close(writerObj);
end