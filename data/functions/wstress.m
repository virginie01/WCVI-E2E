function [taux, tauy, u10, v10] = wstress(u, v, z0)
%WCVIE2E_WSTRESS Compute wind stress using Large & Pond (1981)
%
%   [taux, tauy, u10, v10] = WCVIE2E_wstress(u, v, z0)
%
% DESCRIPTION
%   Computes wind stress and equivalent 10 m wind velocities using the
%   Large and Pond (1981, JPO) bulk aerodynamic formulation.
%
% INPUTS
%   u, v : nt x nx arrays
%       Eastward (u) and northward (v) wind components (m s^-1),
%       measured at height z0.
%
%   z0   : scalar (optional)
%       Height of wind sensor above sea surface (m).
%       Default = 10 m.
%
% OUTPUTS
%   taux, tauy : nt x nx arrays
%       Eastward and northward wind stress components (dyn cm^-2).
%
%   u10, v10   : nt x nx arrays
%       Equivalent wind components at 10 m height (m s^-1).
%
% NOTES
%   - Zero-wind and NaN values are handled safely.
%   - Stress is set to zero where wind speed is zero.
%   - Formulation follows Large & Pond (1981) with bug fix
%     documented by Signell (2005).
%
% REFERENCES
%   Large, W. G., & Pond, S. (1981). Open ocean momentum flux measurements.
%   Journal of Physical Oceanography, 11, 324–336.
%
% Original author: Rich Signell (USGS)
% WCVI-E2E adaptation & cleanup: Virginie Bornarel
%

% -------------------------------------------------------------------------
% Input handling
% -------------------------------------------------------------------------

if nargin < 3 || isempty(z0)
    z0 = 10;  % default measurement height (m)
end

validateattributes(z0, {'numeric'}, {'scalar','positive'});

if ~isequal(size(u), size(v))
    error('WCVIE2E_wstress:DimensionMismatch', ...
        'u and v must have the same dimensions.');
end

[nt, nx] = size(u);

% -------------------------------------------------------------------------
% Physical constants
% -------------------------------------------------------------------------

kappa = 0.41;           % von Kármán constant
rho_a = 1.25e-3;        % air density (g cm^-3)
a1 = (1 / kappa) * log(z0 / 10);

% Preallocate outputs
u10  = nan(nt, nx);
v10  = nan(nt, nx);
taux = nan(nt, nx);
tauy = nan(nt, nx);

% -------------------------------------------------------------------------
% Loop over spatial columns
% -------------------------------------------------------------------------

for ix = 1:nx

    ucol = u(:, ix);
    vcol = v(:, ix);

    % Valid (non-zero, finite) wind entries
    valid = isfinite(ucol) & isfinite(vcol) & ...
            (abs(ucol) + abs(vcol) > 0);

    if ~any(valid)
        continue
    end

    % Complex wind representation
    w = complex(ucol(valid), vcol(valid));
    wspd = abs(w);

    % ---------------------------------------------------------------------
    % Iterative solution for effective wind speed at 10 m
    % ---------------------------------------------------------------------

    d1 = inf(size(wspd));
    c  = zeros(size(wspd));

    while max(abs(c - d1)) > 0.01
        c = d1;

        % Drag coefficient (Large & Pond)
        cd = 1.205e-3 * ones(size(wspd));
        high = c > 11;
        cd(high) = (0.49 + 0.065 .* c(high)) * 1e-3;

        % Adjusted wind speed
        d1 = wspd ./ (1 + a1 .* sqrt(cd));
    end

    % Wind stress magnitude (dyn cm^-2)
    tmag = rho_a .* cd .* d1 * 1e4;

    % 10 m wind components
    w10 = (d1 ./ wspd) .* w;

    % Store results
    u10(valid, ix)  = real(w10);
    v10(valid, ix)  = imag(w10);
    taux(valid, ix) = tmag .* u10(valid, ix);
    tauy(valid, ix) = tmag .* v10(valid, ix);

    % Zero-wind entries explicitly set to zero
    zero = (abs(ucol) + abs(vcol)) == 0;
    u10(zero, ix)  = 0;
    v10(zero, ix)  = 0;
    taux(zero, ix) = 0;
    tauy(zero, ix) = 0;

end

end