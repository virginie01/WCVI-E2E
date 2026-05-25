function [diagnames, A] = setdiagnosticsvars(BioIn, names, ismixed, A)
% SETDIAGNOSTICSVARS Define diagnostic variables for the WCVI-E2E model
%
% INPUTS:
%
%   BioIn   : Structure containing biological input settings and flags
%   names   : Cell array (nbsv x 3) of state variable names (short name, long name, units)
%   ismixed : Logical flag indicating whether the mixed-layer version is used
%   A       : Structure with indices and metadata for state variables (modified in place)
%
% OUTPUTS:
%
% diagnames: Cell array (n x 3) describing all diagnostic variables:
%                 - db/dt for each state variable
%                 - Fluxes between groups
%                 - Predefined diagnostic outputs (e.g. nutrient limitation)
%               Column 1 = short name
%               Column 2 = long description
%               Column 3 = units
%
% A        : Updated structure, now containing:
%              - A.flux : List of all fluxes between groups
%              - A.reroute : User-defined rerouted fluxes (if any)
%              - A.ecosimppflag : Flag for which primary production formulation is used
%
% This file was derived from the original setdiagnosticsvars routine
% developed by Kelly Kearney for the WCE/NEMURO framework and modified
% for the WCVI-E2E coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2014 Kelly Kearney
%
% Modifications and extensions by Virginie Bornarel (2017–2026) include:
%   - revised diagnostic handling for mixed-layer and WCVI-E2E configurations
%   - support for diapause and fisheries-related flux bookkeeping
%   - updated flux-list interfaces and rerouting logic
%   - removal of unused iron-related diagnostics
%   - expanded documentation and diagnostic metadata organization
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.


% --------------------------
% Predefined diagnostic variables (mostly NEMURO)
% --------------------------

diagnames = {...
    'PSpsmax',      'Temp-mediated photsynthesis (small)'    '/s'
    'PLpsmax',      'Temp-mediated photsynthesis (large)'    '/s'
    'PSlightlim',   'Light limitation (small)',         'no units'
    'PLlightlim',   'Light limitation (large)',         'no units'
    'PSno3lim',     'Nitrate limitation (small)',       'no units'
    'PLno3lim',     'Nitrate limitation (large)',       'no units'
    'PSnh4lim',     'Ammonium limitation (small)',      'no units'
    'PLnh4lim',     'Ammonium limitation (large)',      'no units'
    'PLsilim',      'Silica limitation (large)'         'no units'
    'PStemplim'     'Temperature factor (small)',       'no unit'
    'PLtemplim'     'Temperature factor (large)',       'no unit'
    'I',            'Irradiance'                        'W m^-2'
    'kappa',        'Attenuation coefficient'           'm^-1'
    'kp'            'Attenuation self-shading only',    'm^-1'
    'extrasi'       'Extra silica due to negatives'     'mol Si m^-3'};


tf = ~ismember(diagnames(:,1), {'PStemplim','PLtemplim'});

diagnames = diagnames(tf,:);

% --------------------------
% Derivatives: dB/dt for each group
% --------------------------

dbnames = names;
for idb = 1:size(dbnames,1)
    dbnames{idb,1} = ['d' dbnames{idb,1}];
    dbnames{idb,2} = ['dB/dt: ' dbnames{idb,2}];
    dbnames{idb,3} = 'molN m^-3 s^-1';
end

% --------------------------
% Fluxes between groups
% --------------------------
if BioIn.isnem
    
    % Nemuro defaults
    
    fluxlist = listfluxes('nemurokak', ismixed, A.idx, A.diapause);

else
    
    % Wce defaults

    fluxlist = listfluxes('wce', ismixed, A.idx, A.diapause, A.links, A.fisheries);
    
    % Choose the proper primary-production-related fluxes
    
    A.ecosimppflag = BioIn.ecosimpp;
    if A.ecosimppflag
        ismissing = strcmp(fluxlist(:,1), 'gpp') | ...
                    strcmp(fluxlist(:,1), 'exx') | ...
                    strcmp(fluxlist(:,1), 'res');
    else
        ismissing = strcmp(fluxlist(:,1), 'npp');
    end
    fluxlist = fluxlist(~ismissing,:);

end

% --------------------------
% Add user-defined rerouted fluxes
% --------------------------

if isempty(BioIn.reroute)
    A.reroute = cell(0,5);
else
    A.reroute = BioIn.reroute;
    [tf, loc] = ismember(A.reroute(:,2:4), [names(:,1); 'out'; 'fish']);
    if ~all(tf(:))
        error('Unrecognized critter in reroute table');
    end
    A.reroute(:,2:4) = num2cell(loc);
end

fluxlist = [fluxlist; A.reroute(:,[1 2 4])];
nd = size(fluxlist,1);

fluxnames = cell(nd,3);
for id = 1:nd
    fluxnames{id,1} = sprintf('%s_%02d_%02d', fluxlist{id,:});
    fluxnames{id,2} = sprintf('%s: %02d to %02d', fluxlist{id,:});
    fluxnames{id,3} = 'mol m^-3 s^-1';
end

A.flux = fluxlist; % (was nemflux/wceflux in original)

% --------------------------
% Final diagnostic variable list
% --------------------------

diagnames = [dbnames; fluxnames; diagnames];

    
    
