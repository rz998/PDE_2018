function [U] = burgers_upwind_inviscid(L,Ncells,Nsteps)
% This function solves 1D Burgers' equation using the upwind scheme
% under inviscid conditions
%
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


writerObj = VideoWriter('Burgers_upwind_inviscid.avi');
writerObj.FrameRate = 20;

writerObj.Quality = 100;
open(writerObj);

% Main loop
for i = 1:Nsteps
    
    % Find the gradient 
    
    [dU] = minmod(U);
    
    % Compute cell interface values
    UE = U + dU/2;
    
    UW = U - dU/2;
    
    % Compute the squared term
    [dU2] = upwind(U,UE,UW);
    
    % CFL condition (Cmax = 0.1)

    dt = 0.1*dx/max(abs(U));

    % Evolve with respect to time

    U = U - dt*dU2/dx;

    plot(xav,U);

    xlim([0 L]);
    %ylim([4.646e-4 4.658e-4])
    frame = getframe(gcf);
    
    writeVideo(writerObj,frame);


end
close(writerObj);
    
    function [dU2] = upwind(U,UE,UW)
    % This function compute dU^2 with the upwind scheme
    
    dU2 = zeros(1,Ncells);
    
    for k = 1:Ncells
        
        if U(k) > 0
            
            dU2(k) = 0.5*(UE(k)^2-UE(mod((k-2),Ncells)+1)^2);
        
        elseif U(k) < 0
            
            dU2(k) = 0.5*(UW(mod(k,Ncells)+1)^2-UW(k)^2);
           
        else
            
            dU2(k) = 0;
        
        end
     
    end
    
    end  

    function [dU] = minmod(U)
    % This function uses the minmod limiter to find the gradient 
   
    up_dU = zeros(1,Ncells);
    
    down_dU = zeros(1,Ncells);
    
    central_dU = zeros(1,Ncells);
    
    dU = zeros(1,Ncells);
      
    for k = 1:Ncells
        
        up_dU(k) = 2*(U(k)-U(mod(k-2,Ncells)+1));
    
        down_dU(k) = 2*(U(mod(k,Ncells)+1)-U(k));
    
        central_dU(k) = 0.5*((U(mod(k,Ncells)+1))-U(mod(k-2,Ncells)+1));
        
        dU_vector = [up_dU(k),down_dU(k),central_dU(k)];
        
        if min(dU_vector) > 0
            
            dU(k) = min(dU_vector);
        
        elseif max(dU_vector) < 0
            
            dU(k) = max(dU_vector);
        
        else 
            
            dU(k) = 0;
            
        end
     
    end
    
    end    
close(writerObj);
end