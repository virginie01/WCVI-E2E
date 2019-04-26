function [y1, k_1, V1, H1, X1, k_2, V2, H2, X2, k_3, V3, H3, X3, k_4, V4,...
    H4, X4, varargout] = rk4(ode, it, y0, In, Grd, varargin)

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
                                         
opt = varargin(cellfun(@ischar,varargin));

if all(ismember(opt, {'sflux','source'}))
    
[k_1, V1, H1, X1, srf_Tflx1, sol_Tflx1] = ode(x(it), y0, In, Grd, varargin {:});
[k_2, V2, H2, X2, srf_Tflx2, sol_Tflx2] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*k_1), In, Grd, varargin{:});
[k_3, V3, H3, X3, srf_Tflx3, sol_Tflx3] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*k_2), In, Grd, varargin{:});
[k_4, V4, H4, X4, srf_Tflx4, sol_Tflx4] = ode(x(it)+h, y0+(k_3.*h), In, Grd, varargin{:});

varargout{1} = srf_Tflx1;
varargout{2} = sol_Tflx1;
varargout{3} = srf_Tflx2;
varargout{4} = sol_Tflx2;
varargout{5} = srf_Tflx3;
varargout{6} = sol_Tflx3;
varargout{7} = srf_Tflx4;
varargout{8} = sol_Tflx4;

else
    
[k_1, V1, H1, X1] = ode(x(it), y0, In, Grd, varargin {:});
[k_2, V2, H2, X2] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*k_1), In, Grd, varargin{:});
[k_3, V3, H3, X3] = ode(x(it)+(0.5.*h), y0+(0.5.*h.*k_2), In, Grd, varargin{:});
[k_4, V4, H4, X4] = ode(x(it)+h, y0+(k_3.*h), In, Grd, varargin{:});

end

y1  = y0 + (1/6).*(k_1+2.*k_2+2.*k_3+k_4).*h;  % main equation

%% Deal with optional parameter name/value pairs bval and sval

sval.t=Grd.datatime;
sval.x=Grd.boxx;
sval.data=zeros(Grd.datant+1,Grd.nx);

bval.t=Grd.datatime;
bval.x=Grd.boxx;
bval.data=zeros(Grd.datant+1,Grd.nx);

options = struct('sval',sval,'bval',bval);
              
optionsNames=fieldnames(options);

for i=1:length(optionsNames)
if any(strcmp(varargin,optionsNames{i}))
   index = find(strcmp(varargin,optionsNames{i}));
   options.optionsNames{i}=varargin{index+1};           
end
end

if options.sval.data~=zeros(Grd.datant+1,Grd.nx)
   idx=options.sval.t==x(it);
   y1(1,:)=options.sval.data(idx+1,:);
end

if options.bval.data~=zeros(Grd.datant+1,Grd.nx)
    idx=options.bval.t==x(it);
    y1(end,:)=options.bval.data(idx+1,:);
end
