%==========================================================================
% PHYSICAL FORCING VARIABLES — WCVI-E2E MODEL
%==========================================================================
%
% This script prepares the external physical forcing fields required by
% the WCVI-E2E physical–biogeochemical model.
%
% Each forcing variable represents a physically distinct process and is
% provided on either a surface grid, full grid, or boundary interface.
%
% -------------------------------------------------------------------------
% SURFACE FLUXES
% -------------------------------------------------------------------------
% Qi        : Incoming shortwave radiation (W m^-2)
% airtmp    : Near-surface air temperature (°C)
% dewptT    : Dew point temperature (°C)
%
% -------------------------------------------------------------------------
% WIND & UPWELLING
% -------------------------------------------------------------------------
% UWndSpd   : Zonal wind velocity at sensor height (m s^-1)
% VWndSpd   : Meridional wind velocity at sensor height (m s^-1)
% Tauy      : Alongshore wind stress (N m^-2), derived from U/V winds
% Wspd10    : Equivalent neutral-stability wind speed at 10 m (m s^-1)
% Xfil      : Local upwelling index (m s^-1)
% Xfilphy   : Physically filtered upwelling index
%
% -------------------------------------------------------------------------
% VERTICAL MIXING & BUOYANCY
% -------------------------------------------------------------------------
% Mld       : Mixed layer depth (m)
% Entrnmnt : Entrainment rate at the base of the mixed layer (s^-1)
% P         : Precipitation rate (m s^-1)
%
% -------------------------------------------------------------------------
% ALONGSHORE CURRENTS
% -------------------------------------------------------------------------
% CU        : California Undercurrent transport proxy (s^-1)
% DC        : Davidson Current transport proxy (s^-1)
% SBC       : Shelf Break Current transport proxy (s^-1)
%
% -------------------------------------------------------------------------
% LATERAL BOUNDARY CONDITIONS
% -------------------------------------------------------------------------
% tmp       : Boundary temperature forcing (°C)
% sal       : Boundary salinity forcing (psu)
%
% All forcing fields are interpolated to the ODE solver timestep and
% optionally smoothed prior to use in the model.
%
%==========================================================================

%% Directory 

rootdir = fileparts(fileparts(mfilename('fullpath')));
Forcingdir = fullfile(rootdir,'data','physics','Forcing');
LBdir = fullfile(rootdir,'data','physics','LB');

%% Grid parameters to transfer from parseinput.m

In.dz = [31.6349 37.7302; 78.3651 952.2698; 10 10];
In.dx = [40000 40000; 40000 40000; 40000 40000];
In.syear = 1992;
In.eyear = 1992;
In.dt =3*3600;
In.datadt = 0.5* 3600;
In.Lat = 49.5;
In.whgt = 10;

%% Processing physical forcing datasets

In = loadPhysicalForcing(Forcingdir, LBdir);
Grd = buildGrid(In);
Forcing = processPhysicalForcing(In, Grd);