function y = ode4splitsnonneg(odefun,tspan,y0,varargin)
%ODE4SPLITSNONNEG Like ode4, but returns intermediate components, no negative
%
% [y, dy, Splits, Diag, failflag] = ode4splitsnonneg(odefun,tspan,y0,G,P,B, ...)
%
% This function extends the ode4 function to calculate diagnostic variables
% and additive components as the ODE is solved.  For the additive
% components aspect, I assume that the ODE function is of the form dy/dt =
% dy1 + dy2 + dy3 + ....  This is particularly designed for biological
% modules in mixed_layer, where the total change is a sum of processes
% (production, grazing, predation loss, etc).
%
% Input variables:
%
%   odefun:     function handle to ode, of form [db,Splitdb,Diag] =
%               fun(t,b,G,P,B) 
%
%   tspan:      vector of time steps to integrate over
%
%   y0:         3D array of initial conditions
%
%   p#:         additional parameters to pass to the ODE function
%
% Output variables:
%
%   y:          nz x nx x nbsv x 2 array. y0 and new values for next time step 
%
%   dy:         nz x nx x nbsv array. dy over time step in mol N.m-3.s-1
%
%   Splits:     structure containing all types of fluxes in molN.m-3.s-1. Each field is a type of flux 
%               and holds a nbsv+2 x nbsv+2 x nz x nx array that corresponds to the average flux between 
%               a source and a sink group in each box over the current time interval
%
%   Diag:       structure of diagnostic variables at each time step
%
%   failflag:   true if any component becomes negative within the
%               Runge-Kutta calculations.



if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a 3D array of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

f0 = feval(odefun,tspan(1),y0,varargin{:});% mol N.m-3.s-1 nz x nx x nbsv final rate of change 

if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

szy = size(y0); % [nz nx nbsv]
nt = length(tspan);% 2

y = zeros([szy nt]); %[nz nx nbsv 2]
dy = zeros([szy nt-1]); % [nz nx nbsv 1] defined for every interval
y(:,:,:,1) = y0;

for i = 2:nt
    
  ti = tspan(i-1); 
  hi = h(i-1);
  yi = y(:,:,:,i-1);

  f1 = feval(odefun,ti,yi,varargin{:});
  f2 = feval(odefun,ti+0.5*hi,yi+0.5*hi*f1,varargin{:});
  f3 = feval(odefun,ti+0.5*hi,yi+0.5*hi*f2,varargin{:}); 
  f4 = feval(odefun,tspan(i),yi+hi*f3,varargin{:});
 
  dy(:,:,:,i-1) = (1/6)*(f1 + 2*f2 + 2*f3 + f4);%molN.m-3.s-1 nz x nx x nbsv- dy change
  y(:,:,:,i) = yi + dy(:,:,:,i-1).*hi; %molN.m-3 nz x nx x nbsv=new concentration at next time step
  
      
end





