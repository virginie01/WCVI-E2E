function Bio = WCVIE2E_initbiorelax(Grd, Bio, In)
%INITBIORELAX Sets up relaxation for biological variables
%
% Fields added to Bio structure:
% 
%   hasrelax:   nbsv x 1 logical array, true if relaxation data provided
%
%   Relax:      1 x nbsv 1 x 1 structure of data for bio relaxation
%               interpolation.  Fields empty for non-relaxed variables.
%               t:      nt x 1 array, times (Days from simulation start
%                       time) 
%               z:      nz x 1 aray, specifying the box in the z dimension
%               x:      nx x 1 array, specifying the box in the x dimension 
%               data:   nz x nx x nt array, biological relaxation profiles
%               (in units defined bio biological module)
%
%   hasflux:    nbsv x 1 logical array, true if extra flux data provided
%
%   ExtraFlux:  1 x nbsv 1 x 1 structure of data for extra flux
%               interpolation.  Fields empty when biological state 
%               variables don't have extra flux data provided.
%               t:      nt x 1 array, times (Days from simulation start
%                       time) 
%               z:      nz x 1 aray, specifying the box in the z dimension
%               x:      nx x 1 array, specifying the box in the x dimension 
%               data:   nz x nx x nt array, extra flux profiles
%               (in units defined bio biological module)
%
%   

%-------------------------
% Relaxation
%-------------------------

% Determine the names of relaxation input variables (xxxrelax, where xxx
% are the short names of all biological state variables)

relaxvar = cellfun(@(x) [x 'relax'], Bio.names(:,1), 'uni', 0);
inflds = fieldnames(In);

% Indicator variable for relaxation

Bio.hasrelax = ismember(relaxvar, inflds);

% Set up interpolation grids

nb = size(Bio.names, 1);

for ib = 1:nb
    if Bio.hasrelax(ib)
        if isempty(In.(relaxvar{ib}))
            Bio.hasrelax(ib) = false;
        else
        Bio.Relax(ib) = initinterpdata('time and space', In.(relaxvar{ib}), Grd);
        end
    end
end


%-------------------------
% Outside flux
%-------------------------

% Determine the names of relaxation input variables (xxxrelax, where xxx
% are the short names of all biological state variables)

extrafluxvar = cellfun(@(x) [x 'flux'], Bio.names(:,1), 'uni', 0);
inflds = fieldnames(In);

% Indicator variable for relaxation

Bio.hasflux = ismember(extrafluxvar, inflds);

% Set up interpolation grids

for ib = 1:nb
    if Bio.hasflux(ib)
        if isempty(In.(extrafluxvar{ib}))
            Bio.hasflux(ib) = false; 
        else
            Bio.ExtraFlux(ib) = initinterpdata('time and space', In.(extrafluxvar{ib}), Grd);
        end
    end
end
