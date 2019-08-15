function [x,density] = montecarloHS(L,N,n,D,dt,t_final)
%   montecarloHS performs Metropolis Monte Carlo of a 1D model of hard
%   spheres under an external potential of x^2
%  
% 
%   Input arguments
%       L           Length of the box
%       N           Number of hard spheres
%       n           Number of nodes
%       D           Diameter of hard spheres
%       dt          Time interval
%       t_final     Final time
%
%   Output arguments
%       x           Vector of nodes
%       density     Density profile
%

% Set up the nodes 
x = linspace(-L/2,L/2,n);

% Compute the maxmium displacement (in nodes)
dx = dt/(L/n);

% Convert the diameter into nodes
d = round(D/(L/n));
p = randi(n,N,1);

% Initial density profile
while any(pdist(p)<=d)
    p = randi(n,N,1);
end

for i = 1:round(t_final/dt)
    for j = 1:N
    % Attempt to move every sphere
    newp = p;
    newp(j,1) = newp(j,1)+round(dx*randn);
    % Periodic boundary
    if newp(j,1) > n
        newp(j,1) = newp(j,1)-n;
    elseif newp(j,1) < 1
        newp(j,1) = newp(j,1)+n;
    end
    % Compute the potential change
    dU = getPotential(p,newp,x,d,j);
    % Compute the probability of the move
    P = min(1,exp(-dU));
    % Make the move with P
    if P > rand(1)
        p = newp;
    end  
    end
    
  
   
end
density = zeros(n,1);
density(p) = density(p)+1;
end

function [dU] = getPotential(p,newp,x,d,j)
% getPotential calculates the potential energy of the hard spheres

% Compute the external potential
Ue_1 = x(p(j,1))^2;
Ue_2 = x(newp(j,1))^2;
% Compute the internal potential
if any(pdist(p)<=d)
    Ui_1 = Inf;
else 
    Ui_1 = 0;    
end
if any(pdist(newp)<=d)
    Ui_2 = Inf;
else
    Ui_2 = 0;
end
if Ui_1 == Inf && Ui_2 == Inf
    dU = Inf;
else 
    dU = Ui_2+Ue_2-Ui_1-Ue_1;
end
end
