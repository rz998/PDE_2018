function [U] = euler(L,Ncells,Nsteps)
% This function solves 1D Euler equation using the Local Lax-Friedrich
% scheme, where the source term is -\rho x.
%   Input arguments
%       L       The length of interest
%       Ncells  Number of cells
%       Nsteps  Number of time steps
%
%   Output arguments
%       U       Density and momentum profile


% Set the cells
dx = L/Ncells;
x = -L/2:dx:L/2;
xav = 0.5*(x(1:Ncells)+x(2:Ncells+1));

% U is a vector whose first row consists of the density and second row
% consists of the momentum.
U = zeros(2,Ncells);


% Initial condition, p = cos(pi*x/L) and pu = 0.
for j = 1:Ncells

    U(1,j) = integral(@(x) cos(pi.*x/L),x(j),x(j+1))/(dx*integral(@(x) cos(pi.*x/L),-L/2,L/2));

end

writerObj = VideoWriter('euler');
writerObj.FrameRate = 20;

writerObj.Quality = 100;
open(writerObj);


% Main loop
t=0;
for i = 1:Nsteps
    
    % Clear the figure
    clf
    
    % Compute numerical flux and time step
    [dF dt] = LocalLaxFriedrich(U);
    
    % Evolve with respect to time
    U = U + dt*(-dF/dx-[zeros(1,Ncells);U(1,:).*xav]);
    
    t = t+dt;
    
    subplot(2,1,1);
    plot(xav,U(1,:));
    xlim([-L/2 L/2]);
   % ylim([0 1]);
    title('Density profile')
    

    subplot(2,1,2);
    plot(xav,U(2,:));
    xlim([-L/2 L/2]);
   % ylim([-0.6 0.6]);
    title('Momentum profile')
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame)
end


  
   function [dF dt] = LocalLaxFriedrich(U)
   % This function uses the Lax-Friedrich scheme to compute the flux
   
   dF = zeros(2,Ncells);
      
   Fplus = zeros(2,1);
   
   Fminus = Fplus;
    
   for k = 1:Ncells
   
   % Compute the Jacobian matrix
   A = [0,-U(2,k).^2/U(1,k).^2+1,1,2*U(2,k)/U(1,k)];
   
   a = max(abs(A));
   
   dt = 0.1*dx/a;
   
   Fplus(1) = U(2,mod(k,Ncells)+1);
   
   Fplus(2) = U(2,mod(k,Ncells)+1)^2/U(1,mod(k,Ncells)+1) + U(1,mod(k,Ncells)+1);
   
   Fminus(1) = U(2,mod(k-2,Ncells)+1);
   
   Fminus(2) = U(2,mod(k-2,Ncells)+1)^2/U(1,mod(k-2,Ncells)+1) + U(1,mod(k-2,Ncells)+1);
   
   dF(:,k) = 0.5*(Fplus-Fminus)-a*(U(:,mod(k,Ncells)+1)-2*U(:,k)+U(:,mod(k-2,Ncells)+1));
        
   end
    
   end  
close(writerObj);
end