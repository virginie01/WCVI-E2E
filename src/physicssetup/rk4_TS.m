function [y1, dy, Splits] = rk4_TS(ode, it, y0, In, Grd, varargin)

% rk4 calculates ODE using Runge-Kutta 4th order method
%
% INPUTS:
%
% ode:         Function handle to call appropriate ODE.
%
% it:          current time step
%
% y0:          nz x nx array. current conditions for the variable considered.
%
% In:          Structure holding user-supplied input variables.
%
% Grd:         Structure holding spatial and temporal grid data for 
%              WCVIE2E_physicalmodel simulations.
%
% varargin:    additional datasets/parameters needed for ode function. Have to
%              be ordered in the same way as in ode function definition. See
%              ode function definition for description of parameters.
%
% OUTPUTS:  
%
% y1:          nz x nx array. Computed solution. Variable of interest at 
%              next time step.
%
% k_1,...,k_4: Runge Kutta first to fourth order coefficients for diagnostic
%              purposes
%

h = In.dt;                                          % step size (s)
x = Grd.time;                                       % elapsed time (s)
                                         

    
[f1, S1] = ode(x(it), y0, In, Grd, varargin {:});
[f2, S2] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*f1), In, Grd, varargin{:});
[f3, S3] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*f2), In, Grd, varargin{:});
[f4, S4] = ode(x(it)+h, y0+(f3.*h), In, Grd, varargin{:});

dy =(1/6) *(f1 +2.*f2 + 2.*f3 + f4);
y1  = y0 + dy.*h;  


fld = fieldnames(S1);
for is = 1:length(fld)
    Splits.(fld{is}) = (1/6).*(S1.(fld{is})+2.*S2.(fld{is})+2.*S3.(fld{is})+S4.(fld{is}));
end
