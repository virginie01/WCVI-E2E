function [diagnames, A] = setdiagnosticsvars(BioIn, names, ismixed, A)

% OUTPUTS
%
% diagnames: n x 3 cell array describing diagnostic variables from NEMURO, 
%            db/dt for each group, and fluxes in wce model.Column 1 = short 
%            name; column 2 = long name; column 3 = units. 
%
% Diagnostic variables from nemuro

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

% db/dt for each group

dbnames = names;
for idb = 1:size(dbnames,1)
    dbnames{idb,1} = ['d' dbnames{idb,1}];
    dbnames{idb,2} = ['dB/dt: ' dbnames{idb,2}];
    dbnames{idb,3} = 'molN m^-3 s^-1';
end

% Fluxes between groups, by process type, defaults

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

% Add user-rerouted fluxes

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

% Combine all diagnostics variables
% TODO: make consistent between nemurokak and wce?

diagnames = [dbnames; fluxnames; diagnames];

    
    
