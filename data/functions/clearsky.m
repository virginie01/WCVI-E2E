function Q0 = clearsky(start_date, t, Lat)
%WCVIE2E_CLEARSKY Clear-sky shortwave irradiance (Seckel & Beaudry, 1973)
%
%   Q0 = WCVIE2E_clearsky(start_date, t, Lat)
%
% DESCRIPTION
%   Computes clear-sky irradiance following Seckel & Beaudry (1973),
%   as reported in Reed (1977). This formulation is used in the
%   calculation of net longwave heat flux.
%
% INPUTS
%   start_date : 1 x 6 date vector [YYYY MM DD hh mm ss]
%       Reference date corresponding to t = 0.
%
%   t          : vector
%       Time in seconds since start_date.
%
%   Lat        : scalar
%       Latitude in degrees (positive north).
%
% OUTPUT
%   Q0         : vector (same size as t)
%       Clear-sky irradiance (W m^-2).
%
% NOTES
%   - Valid for latitudes between 20°S and 60°N.
%   - Outside this range, coefficients are extrapolated and
%     cloud correction factors may be inaccurate.
%
% REFERENCES
%   Seckel, G. R., & Beaudry, F. H. (1973)
%   Reed, R. K. (1977)
%

% -------------------------------------------------------------------------
% Time of year (days since Jan 1)
% -------------------------------------------------------------------------

% Convert time to fractional day of year (consistent with legacy code)
t_year = datenum(start_date) + t./86400 - ...
         datenum([start_date(1) 1 1 0 0 0]);

% Seasonal phase angle (radians)
phi = pi * ((t_year - 21) * (360/365)) / 180;

% -------------------------------------------------------------------------
% Latitude-dependent coefficients
% -------------------------------------------------------------------------

if (Lat >= 40) && (Lat <= 60)

    A0 = 342.61 - 1.97*Lat - 0.018*(Lat^2);
    A1 = 52.08  - 5.86*Lat + 0.043*(Lat^2);
    B1 = -4.80  + 2.46*Lat - 0.017*(Lat^2);
    A2 = 1.08   - 0.47*Lat + 0.011*(Lat^2);
    B2 = -38.79 + 2.43*Lat - 0.034*(Lat^2);

elseif (Lat >= -20) && (Lat < 40)

    A0 = -15.82 + 326.87*cos(pi*Lat/180);
    A1 = 9.63   + 192.44*cos(pi*(Lat+90)/180);
    B1 = -3.27  + 108.70*sin(pi*Lat/180);
    A2 = -0.64  + 7.80*(sin(pi*(Lat-45)/180))^2;
    B2 = -0.50  + 14.42*(cos(pi*(Lat-5)/180))^2;

else
    warning('WCVIE2E_clearsky:LatitudeOutOfRange', ...
        ['Seckel & Beaudry formulation is only valid between 20S and 60N.\n', ...
         'Results outside this range may lead to incorrect cloud corrections.']);
end

% -------------------------------------------------------------------------
% Clear-sky irradiance
% -------------------------------------------------------------------------

Q0 = A0 + ...
     A1*cos(phi) + B1*sin(phi) + ...
     A2*cos(2*phi) + B2*sin(2*phi);

end