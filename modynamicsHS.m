function [x,density] = modynamicsHS(L,n,N,D,dt,t_final)
%   modynamicsHS performs molecular dynamics of a 1D model of hard
%   spheres under an external potential of x^2
%  
% 
%   Input arguments
%       L           Length of the box
%       n           Number of nodes
%       N           Number of hard spheres
%       D           Diameter of hard spheres
%       dt          Time interval
%       t_final     Final time
%
%   Output arguments
%       density     Density profile
%

% Mass of the hard spheres
Mass = 1;

% Set up the nodes 
x = linspace(-L/2,L/2,n);


% Initialise positions
p = -L/2 + (L/2+L/2)*rand(N,1);
while any(pdist(p)<2^(1/6)*D)
    p = -L/2 + (L/2+L/2)*rand(N,1);
end

% Initialise velocties
v = zeros(N,1);


for i = 1:round(t_final/dt)
    
    p = sort(p);
    
    [F] = getForce(p,N);   
   
    p = p + v*dt + (F/(2*Mass))*dt^2;
    
    [Fnew] = getForce(p,N);

    v = v + ((Fnew + F)/(2*Mass))*dt;
    
    v(abs(v)>abs(p))=0.01*p(abs(v)>abs(p));
    
end

density = zeros(n,1);
pn = round(p*(n-1)/L + (n+1)/2);
density(pn) = density(pn)+1;
plot(x,density);
end

function [F] = getForce(p,N)
% getForce calculates the forces between the hard spheres

% Force due to external potential
F1 = -2*p;

% Force due to internal potential
[dr] = pdist3(p,N);
F2 = (sum(48*dr.^-13-24*dr.^-7))';
F=F1+F2;
end

function [dr] = pdist3(p,N)
% pdist3 calculates the distances between the hard spheres
dr = zeros(N,N);
for i=1:N
    dr(:,i)=p(i)-p;   
end
dr(dr==0) = Inf;
dr(dr>2^(1/6)) = Inf;
dr(dr<-2^(1/6)) = Inf;
end


