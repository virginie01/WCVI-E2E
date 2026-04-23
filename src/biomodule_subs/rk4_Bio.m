function [y1, dy, Splits] = rk4_Bio(odefun, it, y0, In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names)

% rk4 calculates ODE using Runge-Kutta 4th order method
%
% INPUTS:
%
% odefun:      Function handle to call appropriate ODE.
%
% it:          current time step
%
% y0:          nz x nx x nbsv array. Concentration of biological variables in each spatial box at 
%              the current time step (mol N.m-3).
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
% y1:          nz x nx x nbsv array. Computed solution. New concentrations of biological variables 
%              after mixing/advection. Only applicable to groups that are mixed (i.e. "planktonic" groups)
%
% dy:          nz x nx x nbsv array. average dy after mixing/advection occurs over current time 
%              interval in mol N.m-3.s-1
%
% Splits:      structure containing mixing and advection fluxes in molN.m-3.s-1. Each field is a 
%              type of flux and holds a nbsv+2 x nbsv+2 x nz x nx array that corresponds to the 
%              average flux between a source and a sink group in each box over the current time 
%              interval. Here, source/sink groups are the same.
%  
%

h = In.dt;                                          % step size (s)
x = Grd.time;                                       % elapsed time (s)
    
[f1, S1] = feval(odefun, x(it), y0, In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);
[f2, S2] = feval(odefun, x(it)+(0.5.*h), y0+(0.5.*h.*f1), In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);
[f3, S3] = feval(odefun, x(it)+(0.5.*h), y0+(0.5.*h.*f2), In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);
[f4, S4] = feval(odefun, x(it)+h, y0+(f3.*h), In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names);

dy =(1/6) *(f1 +2.*f2 + 2.*f3 + f4);
y1  = y0 + dy.*h; 
%y1 = y0;


fld = fieldnames(S1);
for is = 1:length(fld)
    Splits.(fld{is}) = (1/6).*(S1.(fld{is})+2.*S2.(fld{is})+2.*S3.(fld{is})+S4.(fld{is}));
end

% badthings = isnan(y1) | isinf(y1) | y1<0;

% Check for problems

% if any(badthings(:))
%   [ridx,cidx,tidx] = ind2sub(size(badthings),find(badthings));
%    nb = length(ridx);
%    errstr = cell(nb,1);
%    for ii = 1:nb
%        badbio = y1(ridx(ii), cidx(ii), tidx(ii));
%        if isnan(badbio)
%            errstr{ii} = sprintf('NaN: depth %d, longitude %d, critter %d, time = %d', ridx(ii), cidx(ii), tidx(ii), x(it));
%        elseif isinf(badbio)
%            errstr{ii} = sprintf('Inf: depth %d, longitude %d, critter %d, time = %d', ridx(ii), cidx(ii), tidx(ii), x(it));
%        elseif badbio < 0
%            errstr{ii} = sprintf('Neg: depth %d, longitude %d, critter %d, time = %d', ridx(ii), cidx(ii), tidx(ii), x(it)); 
%        end
%    end
%    errstr = sprintf('  %s\n', errstr{:});
%    errstr = sprintf('Biology out of range:\n%s', errstr);
%    error('rk4_Bio:biologyOutOfRange', errstr);
%end

