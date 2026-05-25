function [y1, dy, Splits] = rk4_TS(ode, it, y0, In, Grd, varargin)
% RK4_TS Runge-Kutta 4th-order integration with diagnostic splits
%
% Computes the next time step of an ODE using the classic Runge-Kutta
% method
%
% INPUTS:
%
% ode:         Function handle to ODE function [dy, splits] = ode(t, y, In, Grd, ...)
%
% it:          Index of current time step
%
% y0:          nz x nx array of state variables at time t.
%
% In:          Structure with simulation input parameters.
%
% Grd:         Structure with temporal/spatial grid data.
%
% varargin:    additional arguments to be passed to the ode function. They 
%              have to be ordered in the same way as in ode function 
%              definition. See ode function definition for description of 
%              parameters.
%
% OUTPUTS:  
%
% y1:          nz x nx array. Updated state variable at next time step (t + dt)
%
% dy:          Rate of change over the time step
%
% Splits:      Struct of diagnotics 
%
%
% Copyright (c) 2026 Virginie Bornarel
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% Time step

h = In.dt;                                          % step size (s)
x = Grd.time;                                       % elapsed time (s)
                                         

% Runge-Kutta steps

[f1, S1] = ode(x(it), y0, In, Grd, varargin {:});
[f2, S2] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*f1), In, Grd, varargin{:});
[f3, S3] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*f2), In, Grd, varargin{:});
[f4, S4] = ode(x(it)+h, y0+(f3.*h), In, Grd, varargin{:});

% Update

dy =(1/6) *(f1 +2.*f2 + 2.*f3 + f4);
y1  = y0 + dy.*h;  

% Compute diagnostic splits

fld = fieldnames(S1);
for is = 1:length(fld)
    Splits.(fld{is}) = (1/6).*(S1.(fld{is})+2.*S2.(fld{is})+2.*S3.(fld{is})+S4.(fld{is}));
end
