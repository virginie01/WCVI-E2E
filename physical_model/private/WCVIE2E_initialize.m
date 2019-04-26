function [Grd, Ht, Wnd, VMxng, TS, LB, P, Arch, Tdiag, Sdiag]=WCVIE2E_initialize(In)
%INITIALIZE Initialize parameters for the physical model
%
% This function sets up most forcing datasets used throughout the WCVIE2E.  
% It also preallocates and sets the initial conditions for many variables 
% that will be modified throughout the model simulation.
%
% Input variables:
%
%   In:     structure holding user-supplied input variables
%
% Output variables:
%
%   Grd:    structure holding temporal and spatial grid parameters:
%
%           z:          nz*nx array, depth coordinate at the center of
%                       each box (negative, m) 
%
%           zp:         (nz+1)*nx array, depth coordinate at the edges of
%                       each box (negative, m) 
%
%           nz:         number of vertical levels
%
%           boxz:       box names in the z dimension (column vector cell
%                       array)
%
%           
%           x:          nz*nx array, distance from the shore at the
%                       center of each box
%
%           xp:         nz*(nx+1) array, distance from the shore at the edges
%                       of each box
%           
%           nx:         number of horizontal levels 
%
%           boxx:       box names in the x dimension (row vector cell
%                       array)
%
%           tmax:       simulation length (seconds)
%
%           nt:         total number of internal time iterations
%
%           time:       1 x nt array, time elapsed from model start time to
%                       the beginning of each time interval (seconds) 
%
%           start_date: 1 x 6 array, date vector for simulation start date.
%                       This will always be Jan 1 of the specified start
%                       year.
%
%           end_date:   1 x 6 array, date vector for simulation end date.
%                       This will always be Dec 31 of the specified end
%                       year.
%
%           datant:     total number of time iterations for forcing
%                       datasets
%
%           datatime:   1 x datant+1 array, time elapsed from model start time 
%                       to the beginning of each time interval in forcing 
%                       datasets (seconds).  
%
%   Ht:     A nested structure holding structures related to heat forcing:
%
%           Qi:        incident solar radiation from Qi_input
%           airtmp:    air temperature from airtmp_input
%           dewptT:    dew point temperature from dewptT_input
%           Qo:        estimate of clear sky iradiance, based on the 
%                      "Smithsonian Formula" from Seckel and Beaudry, as
%                      reported in Reed, 1977, JPO, 7, pp. 482-485. It is
%                      good for latitudes between 20S and 60N.
%
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      datant+1 x 1 vector, time (seconds from simulation 
%                           start time)
%                   x:      a 1 x nx array describing the boxes in the x 
%                           dimension (except for Qo)
%                   data:   datant+1 x nx or datant+1 x 1 array 
%                           (Qo) holding the data
%
%
%   Wnd:    A nested structure holding structures related to wind forcing:
%           Tauy:   S-N wind stress in N.m-2
%           Tauy2:  S-N wind stress in N.m-2 from Ze
%           Wspd10: wind speed at 10m above sea level (m/s)
%   
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      datant+1 x 1 vector, time (seconds from simulation 
%                           start time)
%                   x:      a 1 x nx array describing the boxes in the x 
%                           dimension 
%                   data:   datant+1 x nx array or datant+1 x 1 array (Tauy2) 
%                           holding the data 
%
%   VMxng:  A nested structure holding structures relating to vertical
%           mixing:
%
%           Mld:    Mixed layer depth coming from mld_input
%                       
%                   This structure consists of the following fields:
%                     
%                   t:      datant+1 x 1 vector, time (seconds from simulation 
%                           start time)
%                   x:      a 1 x nx array describing the boxes in the x 
%                           dimension
%                   data:   datant+1 x nx array holding the data
%
%           Entrnmnt: entrainment rate from entrnmnt_input
%
%                     This structure consists of the following fields:
%                     t:    datant+1 x 1 vector, time (seconds from simulation 
%                           start time)
%                     x:    a 1 x nx array describing the boxes in the x 
%                           dimension
%                     data: datant+1 x nx array holding the data
%
%   TS:     structure holding initial values of temperature, salinity, 
%           density and temperature fluxes
%
%           z:          nz x 1 cell array = boxz
%
%           x:          1 x nx cell array = boxx
%
%           T:          nz x nx array, initial temperature profile (deg C)
%
%           S:          nz x nx array, initial salinity profile (psu)
%
%           Sig:        nz x nx array, density profile.  Note that water is
%                       currently treated as incompressible (kg m^-3)
%
%
%   LB:             structure containing boundary variables:
%
%                   tmp:   1 x 1 structure of data holding lateral boundary
%                          temperatures from LBtmp_input
%
%                          t: datant+1 x 1 array, time (seconds from sim start time)
%                          o: 1 x 6 array specifying the origin
%                          data: datant+1 x 6 array holding the data
%
%                   sal:   1 x 1 structure of data holding lateral boundary
%                          salinities from LBsal_input
%
%                          t: datant+1 x 1 array, time (seconds from sim start time)
%                          o: 1 x 6 array specifying the boundary
%                          data: datant+1 x 6 array holding the data
%
%   P:              1 x 1 structure of data holding precipitation fluxes
%                   from p_input
%                  
%                      t: datant+1 x 1 array, time (seconds from sim start time)
%                      x: 1 x nx array describing the boxes in the x dimension
%                      data: datant+1 x nx array holding the data
%   
%   Arch:   structure holding variables related to archive (i.e. output).
%           The archiving period refers to the In.tarch input.
%
%           space:      logical scalar indicating whether archiving outputs
%                       are spatially-resolved.
%
%           startdate:  nbin x 6 array, date vectors corresponding to
%                       start of each archiving time step.
%
%           enddate:    nbin x 6 array, date vectors corresponding to
%                       end of each archiving time step.
%
%           middate:    nbin x 6 array, date vectors corresponding to
%                       middle of each archiving time step.
%
%           fraction:   1 x nt array, fraction that each model time step
%                       contributes to its archiving time step.
%
%           islast:     1 x nt logical array, true if time step
%                       corresponds to the end of an archiving time step.
%
%           bin:        1 x nt array, index of archiving time step to which
%                       each model time step corresponds.
%
%           nbin:       number of archiving time steps
%
%           fileidx:    nt x 2 array, column one holds the index of the
%                       temporary output file to which results will be
%                       written for each time step, column 2 tells the
%                       index of that time step within the file
%
%           isnewfile:  nt x 1 array, true if model time step is the first
%                       to be written to a new temporary output file
% 
%           endtime:    n x 1, end time of each archiving period (seconds)
%
%           filedates:  nfile x 2 array, indices of first and last archive
%                       step included in each temporary output file.
%
%
%  Tdiag and Sdiag: structures used in case diagnostic tests need to be run         
%
% Charlie Stock
% cstock@alum.mit.edu
%
% modified by Kelly Kearney
% modified by Virginie Bornarel

%------------------------------
%Calculate the grid and
%time step properties
%------------------------------

In.dz=In.dz(:,:);
Grd.zp = [0 0; cumsum(In.dz,1)];
Grd.z =(Grd.zp(1:end-1,:)+Grd.zp(2:end,:))./2;
Grd.nz=size(Grd.z,1);
Grd.boxz={'Upper Layer'; 'Lower Layer'; 'Demersal'};

In.dx=In.dx(:,:);
Grd.xp=cumsum(In.dx,2);
Grd.xp=[zeros(size(Grd.xp,1),1),Grd.xp];
Grd.x=(Grd.xp(:,1:end-1)+Grd.xp(:,2:end))./2;
Grd.nx=size(Grd.x,2);
Grd.boxx={'Shelf', 'Slope'};

if isequal(size(In.syear), [1 6])
    Grd.start_date=In.syear;
else
    Grd.start_date=[In.syear 1 1 0 0 0];  % Start at midnight morning Jan 1
end

if isequal(size(In.eyear), [1 6])
    Grd.end_date=In.eyear;
else
    Grd.end_date=[In.eyear 12 31 24 0 0]; 
end

dnstart=datenum(Grd.start_date);
dnend=datenum(Grd.end_date);

Grd.tmax = (dnend - dnstart)*86400; %number of seconds from 1 Jan 1992 at midnight to 1 jan 2017 at midnight
Grd.nt = floor(Grd.tmax/In.dt);

Grd.time=(0:Grd.nt)*In.dt;

% For screen print counter: print progress at the end of each day

daydiff=diff(floor(Grd.time/86400+datenum(Grd.start_date)));
Grd.newday=[true logical(daydiff)];

% time step properties for interpolating forcing datasets for ODE solver

Grd.datant=floor(Grd.tmax./In.datadt);
Grd.datatime=(0:Grd.datant).*In.datadt;

%--------------------------------------------------------------------------
% 1) Interpolate the forcing variables and boundary conditions onto the
%    model grid (i.e.WCVIE2E_initinterpdata)
%--------------------------------------------------------------------------

% Radiation/heat fluxes

Qi = WCVIE2E_initinterpdata('time and surface grid', In.Qi_input, Grd);
airtmp = WCVIE2E_initinterpdata('time and surface grid', In.airtmp_input, Grd);
dewptT = WCVIE2E_initinterpdata('time and surface grid', In.dewptT_input, Grd);

% Wind-related

UWndSpd = WCVIE2E_initinterpdata('time and surface grid', In.uWndSpd_input, Grd);
VWndSpd = WCVIE2E_initinterpdata('time and surface grid', In.vWndSpd_input, Grd);
Tauy2=WCVIE2E_initinterpdata('time', In.tauy2_input, Grd);% data from Ze for comparison purposes

% Vertical mixing-related

Mld = WCVIE2E_initinterpdata('time and surface grid', In.mld_input, Grd);
Entrnmnt = WCVIE2E_initinterpdata('time and surface grid', In.entrnmnt_input, Grd);

% buoyancy fluxes

P= WCVIE2E_initinterpdata('time and surface grid', In.p_input, Grd);

% lateral boundary conditions (i.e. temperature and salinity)

tmp = WCVIE2E_initinterpdata('time and boundaries', In.LBtmp_input, Grd);
sal = WCVIE2E_initinterpdata('time and boundaries', In.LBsal_input, Grd);

%--------------------------------------------------------------------------
% 2) Interpolate forcing time-series data over finer time steps for 
%    ODE solver (i.e. interp1) ---> 30 mins
%--------------------------------------------------------------------------
% Radiation/heat fluxes
Qi.data=interp1(Qi.t,Qi.data,Grd.datatime);
Qi.t=Grd.datatime;

airtmp.data=interp1(airtmp.t,airtmp.data,Grd.datatime);
airtmp.t=Grd.datatime;

dewptT.data=interp1(dewptT.t,dewptT.data,Grd.datatime);
dewptT.t=Grd.datatime;

% Wind-related

UWndSpd.data=interp1(UWndSpd.t,UWndSpd.data,Grd.datatime);
UWndSpd.t=Grd.datatime;

VWndSpd.data=interp1(VWndSpd.t,VWndSpd.data,Grd.datatime);
VWndSpd.t=Grd.datatime;

Tauy2.data=interp1(Tauy2.t,Tauy2.data,Grd.datatime);
Tauy2.t=Grd.datatime;

% Vertical mixing-related

Mld.data=interp1(Mld.t,Mld.data,Grd.datatime);
Mld.t=Grd.datatime;

Entrnmnt.data=interp1(Entrnmnt.t,Entrnmnt.data,Grd.datatime);
Entrnmnt.t=Grd.datatime;

% buoyancy fluxes

P.data=interp1(P.t,P.data,Grd.datatime);
P.t=Grd.datatime;

% lateral boundary conditions (i.e. temperature and salinity)

tmp.data=interp1(tmp.t,tmp.data,Grd.datatime);
tmp.t=Grd.datatime;

sal.data=interp1(sal.t,sal.data,Grd.datatime);
sal.t=Grd.datatime;

%--------------------------------------------------------------------------
% 3) Smooth forcing time-series data ---> moving average over 24hrs-daily
%    mean
%--------------------------------------------------------------------------

span = (24.*3600)./In.datadt;

Tauy2.data=smooth(Tauy2.data,span);

for i=1:Grd.nx
Qi.data(:,i)=smooth(Qi.data(:,i),span);
airtmp.data(:,i)=smooth(airtmp.data(:,i),span);
dewptT.data(:,i)=smooth(dewptT.data(:,i),span);

UWndSpd.data(:,i)=smooth(UWndSpd.data(:,i),span);
VWndSpd.data(:,i)=smooth(VWndSpd.data(:,i),span);

Mld.data(:,i)=smooth(Mld.data(:,i),span);
Entrnmnt.data(:,i)=smooth(Entrnmnt.data(:,i),span);

P.data(:,i)=smooth(P.data(:,i),span);
end

for i=1:6
tmp.data(:,i)=smooth(tmp.data(:,i),span);
sal.data(:,i)=smooth(sal.data(:,i),span);
end
%--------------------------------------------------------------------------
% Some calculations necessary for the heat fluxes
%--------------------------------------------------------------------------
% Clear sky irradiance is based on the "Smithsonian Formula" from Seckel 
% and Beaudry, as reported in Reed, 1977, JPO, 7, pp. 482-485.  It is good 
% for latitudes between 20S and 60N.

Qo.data = WCVIE2E_clearsky(Grd.start_date, Grd.datatime, In.Lat)';
Qo.t = Grd.datatime;
%--------------------------------------------------------------------------
%Calculate surface wind stress parameters
%--------------------------------------------------------------------------
% The wstress routine comes from Rich Signell's RPSstuff toolbox.  
% This gives the stress using the Large and Pond (1981) formula and 
% calculates the 10m wind velocity based on an assumed logarithmic wind 
% velocity profile in a neutral atmosphere.  

%This toolbox can be found at:
%http://woodshole.er.usgs.gov/staffpages/rsignell/rsignell.html 

[~, tauy, uwspeed10, vwspeed10] = WCVIE2E_wstress(UWndSpd.data, VWndSpd.data, In.whgt);

wspeed10 = abs(uwspeed10 + sqrt(-1)*vwspeed10);

% Translate from dynes/cm2 to Newtons/m2

tauy = 0.1*tauy; %datant+1 x nx matrix

% Combine for quicker interpolation later
Tauy.t = UWndSpd.t;
Tauy.x = UWndSpd.x;
Tauy.data = tauy; % ntdata x nx array; column 1 = shelf; column 2 =slope

% Combine for quicker interpolation later
Wspd10.t = UWndSpd.t;
Wspd10.x = UWndSpd.x;
Wspd10.data = wspeed10; % ntdata x nx array; column 1 = shelf; column 2 =slope

%--------------------------------------------------------------------------
%Initialize temperature, salinity, density profiles, and temperature fluxes
%--------------------------------------------------------------------------

TS.z=Grd.boxz;
TS.x=Grd.boxx;

TS.T=In.t_input;

TS.S=In.s_input;

if any(isnan(TS.T(:,:)))
    error('NaN found in initial temperature data');
end
if any(isnan(TS.S(:,:)))
    error('NaN found in initial salinity data');
end

TS.Sig = WCVIE2E_sw_dens0(TS.S, TS.T);

%----------------------------------
% Create structures
%----------------------------------
Ht = struct('Qi',Qi,'airtmp',airtmp,'dewptT',dewptT,'Qo',Qo);

Wnd = struct('Tauy',Tauy,'Tauy2',Tauy2,'Wspd10',Wspd10);

VMxng = struct('Mld',Mld,'Entrnmnt',Entrnmnt);

LB = struct('tmp',tmp,'sal',sal);

%----------------------------------
%Set up archiving
%----------------------------------
% TODO: Three archive variables: beginarchive, tarch, and endarchive
% In parseinput, these need to be set uP
% begin/endarchive = NaN, no cutoff

nout = length(In.tarch);%vector corresponding to number of model ensembles.
                        % Needs to be a vector if I want several output
                        % folders for each model ensemble

%calculate edges of bins for each output file

binedge = cell(nout,1);
for io = 1:nout
    
    %Edges, including entire simulation
    
    binedge{io} = 0:In.tarch(io):max(Grd.time);
    if max (binedge{io}) < max(Grd.time)
        binedge{io} = [binedge{io} max(Grd.time)];
    end
    binedge{io} = binedge{io}(:);
    
     % Adjust for begin/endarchive
    
    if ~isnan(In.beginarchive(io))
        t1 = (In.beginarchive(io) - dnstart)*86400;
        isbefore = binedge{io} < t1;
        binedge{io} = [0; binedge{io}(~isbefore)];  
    end
    if ~isnan(In.endarchive(io))
        t1 = (In.endarchive(io) - dnstart)*86400;
        isafter = binedge{io} > t1;
        binedge{io} = [binedge{io}(~isafter); Grd.tmax];
    end
    
end

% Eliminate any requested outputs that are outside the simulation timespan

isout = cellfun(@(x) length(x)<2, binedge);
if any(isout)
    binedge = binedge(~isout);
    In.outputextension = In.outputextension(~isout);
    nout = length(binedge);
    if nout == 0
        error('None of the specified output file time ranges will result in saved data');
    end
end

% Place timesteps in archiving bins

[n,bin] = deal(cell(nout,1));
for io = 1:nout
    
    [n{io}, bin{io}] = histc(Grd.time, binedge{io});
    
    % Adjust to remove final bin (where timestep equals last bin edge)
    
    n{io}(end-1) = n{io}(end-1) + n{io}(end);
    n{io}(end) = 0;
    bin{io}(bin{io} == length(binedge{io})) = length(binedge{io}) - 1;
end
    
 % Calculate archiving variables

for io = 1:nout
    
    % average over spatial grid
    
    Arch(io).space=In.spatialarch;
    
    % Dates
    
    Arch(io).dateedge = dnstart + binedge{io}./86400;
    Arch(io).middate = (Arch(io).dateedge(1:end-1) + Arch(io).dateedge(2:end))./2;
    
    % For averaging calcs
    
    Arch(io).fraction = 1./n{io}(bin{io});
    Arch(io).islast = [logical(diff(bin{io})) true];
    Arch(io).bin = bin{io};
    Arch(io).nbin = max(Arch(io).bin);
    
	% File names: updated for new folder-based archiving, so .nc extension is 
	% removed.
    
    if nout > 1
        Arch(io).outfile = regexprep(In.outputfile, '.nc', ['_' In.outputextension{io} '.nc']);
    else
        Arch(io).outfile = In.outputfile;
    end
    [pth, fl, ~] = fileparts(Arch(io).outfile);
    Arch(io).outfile = fullfile(pth, fl);
    Arch(io).avg = [];
    Arch(io).ncid = [];
    Arch(io).vid = [];
    
    Arch(io).iens = In.iens;

end
   
%--------------------------------------------------------------------------
% OPTIONAL: create folder and empty files if diagnostics need to be run
%--------------------------------------------------------------------------
    
    Tdiag.folder='T_diagnostic';
    
    Tdiag.file=[];
    Tdiag.ncid=[];
    Tdiag.vid=[];
    
    Sdiag.folder='Sal_diagnostic';
    
    Sdiag.file=[];
    Sdiag.ncid=[];
    Sdiag.vid=[];
    

      



