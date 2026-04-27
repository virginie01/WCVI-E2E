%==========================================================================
% BIOLOGICAL FORCING VARIABLES — WCVI-E2E MODEL
%==========================================================================
%
% This script prepares the external biological forcing fields required by
% the WCVI-E2E coupled physical–biogeochemical–trophic model.
%
% These forcings are used when biological state variables are prescribed
% externally or when boundary conditions are required for planktonic
% tracers and higher trophic levels.
%
% -------------------------------------------------------------------------
% DISSOLVED OXYGEN
% -------------------------------------------------------------------------
% O2
%   Dissolved oxygen concentration (mol O2 m^-3) prescribed over the full
%   spatial grid (x, z) and time.
%
%   This forcing is used when oxygen is not explicitly simulated by the
%   biogeochemical module and is required for processes such as:
%     - ammonification
%     - nitrification
%     - denitrification
%
% -------------------------------------------------------------------------
% FISHERIES TIME SERIES
% -------------------------------------------------------------------------
% BIOTS
%   Fisheries extraction time series following the Ecosim convention.
%
%   Includes catch or fishing mortality rates for selected functional
%   groups, provided at the physical model timestep.
%
%   Each time series is associated with:
%     - a functional group
%     - a unit (e.g. t yr^-1 or t km^-2)
%     - a pool and extraction type
%
% -------------------------------------------------------------------------
% BIOLOGICAL LATERAL BOUNDARY CONDITIONS
% -------------------------------------------------------------------------
% LB.<name>
%   Time-varying lateral boundary concentrations for planktonic and
%   dissolved biogeochemical tracers (mol N m^-3).
%
%   These forcings prescribe concentrations entering the model domain from:
%     - open-ocean boundaries
%     - riverine inputs
%     - alongshore currents (CU, DC, SBC)
%
%   Boundary forcings are applied only to biologically mixed tracers and
%   are repeated climatologically when provided for a single year.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% - All biological forcing datasets are interpolated to the ODE solver
%   timestep prior to model integration.
%
% - Physical meanings, units, and data provenance for each forcing
%   variable are documented in the corresponding data preparation scripts.
%
%==========================================================================


%% Directory 

rootdir = fileparts(fileparts(mfilename('fullpath')));
Forcingdir = fullfile(rootdir,'data','bio','Forcing');
LBdir = fullfile(rootdir,'data','bio','LB');

%% Grid parameters to transfer from parseinput.m

In.dz = [31.6349 37.7302; 78.3651 952.2698; 10 10];
In.dx = [40000 40000; 40000 40000; 40000 40000];
In.syear = 1992;
In.eyear = 1992;
In.dt =3*3600;
In.datadt = 0.5* 3600;

%% Processing biological forcing datasets
Grd = buildGrid(In);
BioIn = loadBiologicalForcing(Forcingdir, LBdir);

% -------------------------------------------------------------------------
% Oxygen forcing (3D, full grid)
% -------------------------------------------------------------------------
BioForcing.O2 = buildO2Forcing(BioIn.o2_input, Grd);

% -------------------------------------------------------------------------
% Fisheries time series (Ecosim-style)
% -------------------------------------------------------------------------
BioForcing.BIOTS = buildFisheriesTS(BioIn.fisheriesTS_input, Grd);

% -------------------------------------------------------------------------
% Lateral boundary conditions for biological tracers
% -------------------------------------------------------------------------
BioForcing.LB = buildBiologicalLateralBC(BioIn.LB, Grd);


