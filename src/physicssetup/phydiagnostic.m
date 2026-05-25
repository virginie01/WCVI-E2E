function [diagnames,flux] = phydiagnostic()
% PHYDIAGNOSTIC  Returns standard physical diagnostics and flux link definitions.
%
%   [diagnames, flux] = PHYDIAGNOSTIC() generates two outputs:
%
%   diagnames:
%       A cell array of physical diagnostic variable names and units.
%       Each row is {name, description, unit}.
%       - Includes state variable rates of change (e.g., temperature, salinity)
%         with units in deg C s^-1 or psu s^-1.
%       - Includes standard intermediate fluxes (e.g., cross-shore advection, heat flux, etc.)
%         formatted as '<fluxname>_<from>_<to>' with units in deg C s^-1 or psu s^-1
%
%   flux:
%       A cell array of flux link definitions.
%       Each row is {fluxname, source variable, sink variable}.
%
% Copyright (c) 2026 Virginie Bornarel
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.


%% --- Define state diagnostic names ---
dbnames = {...
 'dtemp'     'dB/dt:temperature'       'deg C s^-1'
 'dsal'      'dB/dt:salinity'          'psu s^-1'
 };

%% --- Define flux linkages ---

        A.V = [...
            1     1
            2     2];
        
        A.H = [...
            1     1
            2     2];

        
        A.X = [...
            1     1
            2     2];

        
        A.CS = [...
            1     1
            2     2];

        
        A.P = [...
            1     1
            2     2];

        
        A.R = [...
            1     1
            2     2];


        A.VICC = [...
            1     1
            2     2];

        
        A.DC = [...
            1     1
            2     2];


        A.SBC = [...
            1     1
            2     2];

        
        A.CU = [...
            1     1
            2     2];

        
        A.sol_Tflx = [...
            1     1];

        
        A.srf_Tflx = [...
            1     1];
        
        flux = cell(0,3);

%% --- Build flux cell array ---
fld = fieldnames(A);
for ifld = 1:length(fld)
    nlnk = size(A.(fld{ifld}),1);
    
    flux = [flux; repmat(fld(ifld), nlnk, 1), num2cell(A.(fld{ifld}))];
end

nd = size(flux,1);

%% --- Format flux names ---
fluxnames = cell(nd,3);
for id = 1:nd
    fluxnames{id,1} = sprintf('%s_%02d_%02d', flux{id,:});
    fluxnames{id,2} = sprintf('%s: %02d to %02d', flux{id,:});
    fluxnames{id,3} = 'deg C/psu s^-1';
end

%% --- Combine all diagnostics ---
diagnames = [dbnames; fluxnames];


end

