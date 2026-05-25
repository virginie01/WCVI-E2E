function [y1, dy, Splits] = rk4_Bio(odefun, it, y0, In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names)
% RK4_BIO Integrates a biological model using the 4th-order Runge-Kutta method
%
% This function performs one Runge-Kutta 4th-order integration step for a 
% biological tracer model defined on a spatial grid (e.g., WCVIE2E model).
%
% INPUTS:
%   odefun   - Function handle to the ODE solver, must return:
%                dy: rate of change (same size as y0)
%                Splits: structure of flux terms (optional)
%
%   it       - Current time step index (scalar)
%   y0       - nz x nx x nbsv array of biological tracer concentrations at time step it (mol N m⁻³)
%   In       - Structure with user-defined simulation parameters (e.g., In.dt)
%   Grd      - Structure with grid and time data (Grd.time is a 1 x nt array)
%   ismixed  - Logical mask or flag indicating which biological groups are mixed
%   mld      - Mixed layer depth (may vary over space/time)
%   ent, lb, Sig, tauy, xfil, buoy, names
%            - Additional physical or biological forcing variables (specific to model)
%
% OUTPUTS:
%   y1       - nz x nx x nbsv array. Updated concentrations at t+dt
%   dy       - nz x nx x nbsv array. Average rate of change during [t, t+dt] (mol N m⁻³ s⁻¹)
%   Splits   - Struct containing averaged mixing/advection fluxes (mol N m⁻³ s⁻¹).
%              Each field is a (nbsv+2) x (nbsv+2) x nz x nx array representing source-to-sink fluxes.
%
%  
% Copyright (c) 2026 Virginie Bornarel
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.


h = In.dt;                                          % step size (s)
x = Grd.time;                                       % elapsed time (s)
    
[f1, S1] = feval(odefun, x(it), y0, In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);
[f2, S2] = feval(odefun, x(it)+(0.5.*h), y0+(0.5.*h.*f1), In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);
[f3, S3] = feval(odefun, x(it)+(0.5.*h), y0+(0.5.*h.*f2), In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);
[f4, S4] = feval(odefun, x(it)+h, y0+(f3.*h), In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);

dy =(1/6) *(f1 +2.*f2 + 2.*f3 + f4);
y1  = y0 + dy.*h; 


fld = fieldnames(S1);
for is = 1:length(fld)
    Splits.(fld{is}) = (1/6).*(S1.(fld{is})+2.*S2.(fld{is})+2.*S3.(fld{is})+S4.(fld{is}));
end
