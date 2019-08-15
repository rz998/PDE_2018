function [U] = burgers_LF(L,Ncells,Nsteps)
% This function solves 1D Burgers' equation using the Lax-Friedrich scheme
% under inviscid condtions.
%   Input arguments
%       L       The length of interest
%       Ncells  Number of cells
%       Nsteps  Number of time steps
%
%   Output arguments
%       U       Velocity profile

% Set the cells
dx = L/Ncells;
x = 0:dx:L;
xav = 0.5*(x(1:Ncells)+x(2:Ncells+1));
U = zeros(1,Ncells);

% Initial condition
for j = 1:Ncells

    U(j) = integral(@(x) 4+2*cos(2*pi.*x),x(j),x(j+1))/dx;

end

writerObj = VideoWriter('Burgers_LF');
writerObj.FrameRate = 20;

writerObj.Quality = 100;
open(writerObj);

% Main loop
for i = 1:Nsteps
    
    % Compute the numerical flux and time step
    
    [dF dt] = LocalLaxFriedrich(U);
    
    % Evolve with respect to time

    U = U - dt.*dF/dx;

    plot(xav,U);

    xlim([0 L]);
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    
end

    function [dF dt] = LocalLaxFriedrich(U)
    % This function uses the Lax-Friedrich scheme to compute the flux
     
    % Numerical viscosity
    a = abs(U);
    
    dt = 0.5*dx./a;
    
    dF = zeros(1,Ncells);
    
    for k = 1:Ncells
        
        dF(k) =  0.5*(0.5*(U(mod(k,Ncells)+1)^2-U(mod(k-2,Ncells)+1)^2))...
                 +a(k)*(2*U(k)-U(mod(k,Ncells)+1)-U(mod(k-2,Ncells)+1)); 
    
    end
    
    end    
close(writerObj);
end