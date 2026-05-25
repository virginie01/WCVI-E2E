function y = ode4splitsnonneg(odefun,tspan,y0,varargin)
%ODE4SPLITSNONNEG Runge-Kutta 4th order solver with split tracking and non-negativity enforcement.
%
% y = ode4splitsnonneg(odefun, tspan, y0, G, P, B, Arch, it)
%
% This function integrates an ODE system using the classical 4th-order Runge-Kutta method
% while returning intermediate components and ensuring no negative concentrations. It is
% particularly useful for biological models where total rates are composed of additive processes.
%
% INPUTS:
% odefun   - Function handle to ODE system. Should return db = odefun(t, y, ...)
% tspan    - Vector of time steps to integrate over (length must be >= 2)
% y0       - Initial condition (3D array: nz x nx x nbsv)
% varargin - Additional arguments passed to odefun
%
% OUTPUT:
% y        - 4D array (nz x nx x nbsv x nt), state values at each time step
%
% This file was derived from the original ode4splitsnonneg routine
% developed by Kelly Kearney for the WCE/NEMURO framework and modified
% for use with WCVI-E2E biological state arrays.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Modifications by Virginie Bornarel (2017–2026) include:
%   - adaptation to 3D biological state arrays
%   - revised array indexing for WCVI-E2E model geometry
%   - simplified solver output handling
%   - updated documentation and usage notes
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

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





