function [Qi, Qs, Ql, Qlw] = WCVIE2E_calcheat(Qi, atmp, dpnt, Tsrf, wspeed10, Qo, meanQi, alb, nx)
%CALCHEAT Calculate heating factors
%
% [Qi, Qs, Ql, Qlw] = calcheat(Qi, atmp, dpnt, Tsrf, wspeed10, ...
%                                       Qo, meanQi, alb)
%
% Input variables:
%
%   Qi:         1 x nx array. Incident solar radiation (W m^-2). 
%
%   atmp:       1 x nx array. Air temperature (deg C). 
%
%   dpnt:       1 x nx array. Dewpoint temperature (deg C). 
%
%   Tsrf:       1 x nx array. Sea surface temperature (deg C). 
%
%   wspeed10:   1 x nx array. Wind speed 10 m above ocean surface (m/s). 
%
%   Qo:         Scalar. Clear-sky irradiance (W m^-2).
%
%   meanQi:     1 x nx array. Mean observed daily irradiance (W m^-2).
%
%   alb:        scalar. Albedo, i.e. fraction of incoming radiation
%               reflected from the sea surface.
%
%   nx:         Scalar. Number of horizontal levels. 
%
% Output variables:
%
%   Qi:         1 x nx array. Incident solar radiation (W m^-2).
%
%   Qs:         1 x nx array. Sensible heat flux (W m^-2).
%
%   Ql:         1 x nx array. Latent heat flux (W m^-2).
%
%   Qlw:        1 x nx array. Longwave heat flux (W m^-2).
%
%-----------------------------
% Calculate sensible heat flux
% (Qs) based on empirical 
% relationships of Freihe and 
% Schmitt, 1976, JPO, 6, pp. 
% 801-809
%-----------------------------

rho_air = 1.2;  % approximate density of dry air value, kg/m3 (20 C). 
Cp_air = 1000;  % approximate specific heat of dry air J kg-1 K-1
delt = (Tsrf - atmp);
UT = wspeed10.*delt;

% Use relationships on p. 807 of Friehe and Schmitt to calculate the
% vertical velocity, temperature covariance.  

for i = 1:nx
    if UT(i) <= 0
       WT(1,i) = 0.0026 + 0.86e-3*UT(i);
    elseif ((UT(i) > 0) & (UT(i) < 25))   % unstable boundary layers
       WT(1,i) = 0.002 + 0.97e-3*UT(i);
    else                                  % highly unstable boundary layers
       WT(1,i) = 1.46e-3*UT(i);
    end
end

% Calculate Qs w/relationship (1) of F&S, flux is positive if heat is
% moving from the atmosphere to the ocean.

Qs = (-rho_air*Cp_air).*WT;

%-----------------------------
% Calculate latent heat flux 
% (Ql) based on empirical 
% relationships of Freihe and 
% Schmitt, 1976, JPO, 6, pp. 
% 801-809     
%-----------------------------

% Calculate the vapor pressure at each dew point using equation 1-4 of
% Marshall and Plumb (2008, p. 6).  This is an approximation of the
% Clausius-Clapeyron relationship at typical atmospheric temperatures.
A = 611;        % Pascals
Beta = 0.067;   % Celsius^-1

wvp_air = A.*exp(Beta.*dpnt);     % water vapor pressure in air (Pascals)
wvp_srf = A.*exp(Beta.*Tsrf);     % water vapor pressure at water surface
                                  % (Pascals)

R = 8.314;          % universal gas constant J K-1 mol-1
Le = 2.25e3;        % Latent heat of evaporation, J g-1, Marshall and
                    % Plumb, p. 166
                    
dpntK = dpnt+273.15;    % dew point temperature in Kelvin
wtmpK = Tsrf+273.15;    % Air temperature in Kelvin

% Calculate the number grams of water vapor per cubic meter using the 
% ideal gas law, PV = nRT, solve for number of moles (n) and multiply 
% by 18g of water vapor per mole.

wv_air = 18.*wvp_air./(R.*dpntK);      % water vapor, g/m^3
wv_srf = 18.*wvp_srf./(R.*wtmpK);      % water vapor, g/m^3

% calculate the covariance between the vertical velocity and the water
% vapor density (WQ) using Friehe and Schmitt relationship (p. 808)
% (m s-1) (g m-3)

WQ = 1.32e-3.*wspeed10.*(wv_srf-wv_air);

% Multiply WQ by the latent heat of evaporation (J g-1) to get the latent
% heat flux (J s-1 m-2 = watts m-2), the negative gets flux from the
% surface to the atmosphere.

Ql = -Le.*WQ;

%-----------------------------
% Calculate the net longwave 
% heatflux (Qlw) w/the Efimova 
% formula as reported by 
% Simpson and Paulson, 1979, 
% Quart. J. R. Meteor. Soc. 
% 105, pp. 487-502.  
%-----------------------------

% Qlw_clear is the estimated mean daily clear sky irradiance.  meanQi is
% the mean observed irradiance.  The ratio of these terms provides a
% correction for the Qlw due to clouds. 

emis = 0.97;         % emissivity
ea = wvp_air./10^2;   % atmospheric vapor pressure in millibars: 1 bar =
                     % 10^5 pascals, 1 millibar = 100 pascals
boltz = 5.67e-8;     % Stefan-Boltzmann constant [W m2 K-4] 
Qlw_clear  = (emis*boltz).*(wtmpK.^4).*(0.254-0.00495.*ea);

% Find the average irradiance in the surrounding 24 hours

cldfac = max(1-(meanQi./Qo),[0 0]);

% Simple cloud correction factor of Simpson and Paulson

Qlw = -Qlw_clear.*(1 - 0.8.*cldfac);



