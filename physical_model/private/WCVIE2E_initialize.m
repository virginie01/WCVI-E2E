function [Grd, Ht, Wnd, VMxng, TS, LB, Arch]=WCVIE2E_initialize(In)
%INITIALIZE Initialize parameters for the physical model
%
% This function initializes the values of most constant parameters used
% throughout the mixed-layer model.  It also preallocates and sets the
% initial conditions for many variables that will be modified throughout
% the model simulation.
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
%                       each box (m) 
%
%           zp:         (nz+1)*nx array, depth coordinate at the edges of
%                       each box (m) 
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
%           xp:         nz*nx array, distance from the shore at the edges
%                       of each box
%           
%           nx:         number of horizontal levels 
%
%           boxx:       box names in the x dimension (row vector cell
%                       array)
%
%           tmax:       simulation length (days)
%
%           nt:         total number of internal time iterations
%
%           time:       1 x nt array, time elapsed from model start time to
%                       the beginning of each time interval (days) 
%
%           start_date: 1 x 3 array, date vector for simulation start date.
%                       This will always be Jan 1 of the specified start
%                       year.
%
%           end_date:   1 x 3 array, date vector for simulation end date.
%                       This will always be Dec 31 of the specified end
%                       year.
%
%   Ht:     A nested structure holding structures related to heat forcing:
%
%           Qi:     incident solar radiation from Qi_input
%           Airtmp: Air temperature from airtmp_input
%           DewptT: Dew point temperatures from dewptT_input
%           Qo:     Estimate of clear sky irradiance, based on the 
%                   "Smithsonian Formula" from Seckel and Beaudry, as 
%                   reported in Reed, 1977, JPO, 7, pp.482-485.  It is good 
%                   for latitudes between 20S and 60N.
%           MeanQi: nt x 1 array, mean observed daily irradiance (W.m-2)
%
%
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      nt x 1 vector, time (days from simulation start 
%                           time)
%                   x:      a 1 x nx array describing the boxes in the x 
%                           dimension (except for Q0)
%                   data:   nt x nx or nt x 1 (Qo) array holding the data
%
%
%   Wnd:    A nested structure holding structures related to wind forcing:
%
%           UWndSpd:E-W wind speed from uWndSpd_input
%           VWndSpd:N-S wind speed from vWndSpd_input
%           WndSpd10: speed velocity (m.d-1) 10 m above sea level
%   
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      nt x 1 vector, time (days from simulation start 
%                           time)
%                   x:      a 1 x nx array describing the boxes in the x 
%                           dimension 
%                   data:   nt x nx array holding the data
%
%           Wndl:   Local longshore wind stress from wndl_input
%           Wndr:   Remote wind forcing from wndr_input
%
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      nt x 1 vector, time (days from simulation start 
%                           time)
%
%                   data:   nt x 1 array holding the data
%
%   VMxng:  A nested structure holding structures relating to vertical
%           mixing:
%
%           Mld:    Mixed layer depth coming from mld_input
%                       
%                   This structure consists of the following fields:
%                     
%                   t:      nt x 1 vector, time (days from simulation start 
%                           time)
%                   x:      a 1 x nx array describing the boxes in the x 
%                           dimension
%                   data:   nt x nx array holding the data
%
%           Entrnmnt: entrainment rate from entrnmnt_input
%
%                     This structure consists of the following fields:
%                     t:    nt x 1 vector, time (days from simulation start 
%                           time)
%                     x:    a 1 x nx array describing the boxes in the x 
%                           dimension
%                     data: nt x nx array holding the data
%
%   TS:     structure holding variables related to temperature and salinity
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
%           Srelax:     1 x 1 structure of data for salt relaxation
%                       interpolation.  Not included if no salt relaxation
%                       data was provided.
%                       t:    nt x 1 vector, time (days from sim start time) 
%                       z:    nz x 1 array describing the boxes in the z
%                             dimension
%                       x:    1 x nx array describing the boxes in the x 
%                             dimension
%               
%                       data: nz x nx x nt array, salt relaxation profiles
%                             (psu) 
%
%           Trelax:     1 x 1 structure of data for temperature relaxation
%                       interpolation.  Not included if no temperature
%                       relaxation data was provided.
%                       t:    nt x 1 vector, time (days from sim start time) 
%                       z:    nz x 1 array describing the boxes in the z
%                             dimension
%                       x:    1 x nx array describing the boxes in the x 
%                             dimension
%               
%                       data: nz x nx x nt array, temperature relaxation
%                             profiles(deg C) 
%
%
%   LB:             structure containing boundary 
%                   variables:
%                   
%                   physics: structure contining the following fields:
%
%
%                   LBtmp: 1 x 1 structure of data holding lateral boundary
%                          temperatures from LBtmp_input
%
%                          t: nt x 1 array, time (days from sim start time)
%                          o: 1 x 5 array specifying the boundary
%                          data: nt x 5 array holding the data
%
%                   LBsal: 1 x 1 structure of data holding lateral boundary
%                          salinities from LBsal_input
%
%                          t: nt x 1 array, time (days from sim start time)
%                          o: 1 x 5 array specifying the boundary
%                          data: nt x 5 array holding the data
%
%   
%   Arch:   structure holding variables related to archive (i.e. output).
%           The archiving period refers to the In.tarch input.
%
%           startdate:  nbin x 3 array, date vectors corresponding to
%                       start of each archiving time step
%
%           enddate:    nbin x 3 array, date vectors corresponding to
%                       end of each archiving time step
%
%           middate:    nbin x 3 array, date vectors corresponding to
%                       middle of each archiving time step
%
%           fraction:   1 x nt array, fraction that each model time step
%                       contributes to its archiving time step
%
%           islast:     1 x nt logical array, true if time step
%                       corresponds to the end of an archiving time step
%
%           bin:        1 x nt array, index of archiving time step to which
%                       each model time step corresponds
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
%           endtime:    n x 1, end time of each archiving period (s)
%
%           filedates:  nfile x 2 array, indices of first and last archive
%                       step included in each temporary output file
%           
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
Grd.zp = -[0 0; cumsum(In.dz,1)];
Grd.z =(Grd.zp(1:end-1,:)+Grd.zp(2:end,:))./2;
Grd.nz=size(Grd.z,1);
Grd.boxz={'Upper Layer'; 'Lower Layer'; 'Demersal'};

In.dx=In.dx(:,:);
Grd.xp=cumsum(In.dx,2);
Grd.xp=[zeros(size(Grd.xp,1),1),Grd.xp];
Grd.x=(Grd.xp(:,1:end-1)+Grd.xp(:,2:end))./2;
Grd.nx=size(Grd.x,2);
Grd.boxx={'Shelf', 'Slope'};

if isequal(size(In.syear), [1 3])
    Grd.start_date=In.syear;
else
    Grd.start_date=[In.syear,1,1];
end

if isequal(size(In.eyear), [1 3])
    Grd.end_date=In.eyear;
else
    Grd.end_date=[In.eyear,12,31]; 
end

dnstart=datenum(Grd.start_date);
dnend=datenum(Grd.end_date);

Grd.tmax = (dnend - dnstart);
Grd.nt = floor(Grd.tmax/In.dt);

Grd.time=(0:Grd.nt)*In.dt;

% For screen print counter: print progress at the end of each day

daydiff=diff(floor(Grd.time+datenum(Grd.start_date)));
Grd.newday=[true logical(daydiff)];

% Hotstart save index

if ~isempty(In.hotstartdn)
    Grd.savehot=find(Grd.time<(In.hotstartdn-dnstart),1,'last');
else
    Grd.savehot=NaN;
end

%----------------------------------
% Interpolate the forcing variables 
% onto the model grid
%-----------------------------------

Qi = WCVIE2E_initinterpdata('time and surface grid', In.Qi_input, Grd);
Airtmp = WCVIE2E_initinterpdata('time and surface grid', In.airtmp_input, Grd);
DewptT = WCVIE2E_initinterpdata('time and surface grid', In.dewptT_input, Grd);

UWndSpd = WCVIE2E_initinterpdata('time and surface grid', In.uWndSpd_input, Grd);
VWndSpd = WCVIE2E_initinterpdata('time and surface grid', In.vWndSpd_input, Grd);
Wndl = WCVIE2E_initinterpdata('time', In.wndl_input, Grd);
Wndr = WCVIE2E_initinterpdata('time', In.wndr_input, Grd);

Mld = WCVIE2E_initinterpdata('time and surface grid', In.mld_input, Grd);
Entrnmnt = WCVIE2E_initinterpdata('time and surface grid', In.entrnmnt_input, Grd);

VMxng = struct('Mld',Mld,'Entrnmnt',Entrnmnt);

%-------------------------------------
%Initialize temperature, salinity and 
%density profiles. 
%-------------------------------------

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

%--------------------------------------
%set lateral boundary conditions for
%temperature and salinity
%--------------------------------------

LBtmp = WCVIE2E_initinterpdata('time and boundaries', In.LBtmp_input, Grd);
LBsal = WCVIE2E_initinterpdata('time and boundaries', In.LBsal_input, Grd);

LB.physics=struct('LBtmp',LBtmp,'LBsal',LBsal);


%-----------------------------------
%Set up temperature and 
%salinity relaxation
%-----------------------------------

  %salinity
  if In.hassrelax  
      TS.Srelax = WCVIE2E_initinterpdata('time and full grid', In.srelax, Grd);
  end
  %temperature
  if In.hastrelax
      TS.Trelax = WCVIE2E_initinterpdata('time and full grid', In.trelax, Grd);
  end
 
%-----------------------------------
%Calculate surface wind
%stress parameters
%-----------------------------------
% Calculate wind speed at 10m (Wnd.wspeed10).  The quantities are all used 
% for the heat flux calculations.  The wstress routine comes from Rich 
% Signell's RPSstuff toolbox.  This gives the stress using the Large and 
% Pond (1981) formula and calculates the 10m wind velocity based on an 
% assumed logarithmic wind velocity profile in a neutral atmosphere.  

%This toolbox can be found at:
%http://woodshole.er.usgs.gov/staffpages/rsignell/rsignell.html 

[uwspeed10, vwspeed10] = WCVIE2E_wstress(UWndSpd.data, VWndSpd.data, In.whgt);

wspeed10 = abs(uwspeed10 + sqrt(-1).* vwspeed10);%wspeed10 nt x nx matrix

% Combine for quicker interpolation later
WndSpd10.t = UWndSpd.t;
WndSpd10.x = UWndSpd.x;
WndSpd10.data = wspeed10; % nt x nx model; column 1 = shelf; column 2 =slope

%Creation of a nested structure containing all variables related to wind

Wnd = struct('UWndSpd',UWndSpd,'VWndSpd',VWndSpd,'WndSpd10',WndSpd10,...
    'Wndl',Wndl,'Wndr',Wndr);

%--------------------------
% Some calculations 
% necessary for the heat 
% fluxes                
%--------------------------

% Construct a time series of the daily mean irradiance.  This quantity will
% be used in conjunction with Ht.Qo to calculate a cloud correction factor if 
% hswitch = 2.  It is also used in the biological calculations produce PAR
% averaged over 14 daylight hours.

% If you are deriving Qs, Ql, and Qlw from other meteorological inputs an
% estimate of the clear sky irradiance will be needed.  This is based on
% the "Smithsonian Formula" from Seckel and Beaudry, as reported in Reed, 
% 1977, JPO, 7, pp. 482-485.  It is good for latitudes between 20S and 60N.

Qo.t=Qi.t;
Qo.data = WCVIE2E_clearsky(Grd.start_date,Grd.time,In.Lat)';

% calculate the mean observed daily irradiance (watts/m2).  This will
% be used in conjunction with Ht.Qo to calculate a cloud correction
% factor if hswitch = 2.  It is also used in the biological
% calculations.

qi=interp1(Qi.t, Qi.data, Grd.time);% Qi.data = nt x nx array gives qi with same size

MeanQi.t = Qi.t;
MeanQi.x = Qi.x;
for i=1:Grd.nx
MeanQi.data(:,i) = smooth(Grd.time', qi(:,i), 8, 'moving');
end

Ht = struct('Qi',Qi,'Airtmp',Airtmp,'DewptT',DewptT,...
    'Qo',Qo,'MeanQi',MeanQi);
  
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
        t1 = (In.beginarchive(io) - dnstart);
        isbefore = binedge{io} < t1;
        binedge{io} = [0; binedge{io}(~isbefore)];  
    end
    if ~isnan(In.endarchive(io))
        t1 = (In.endarchive(io) - dnstart);
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
    
    % Dates
    
    Arch(io).dateedge = dnstart + binedge{io};
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
    [pth, fl, ex] = fileparts(Arch(io).outfile);
    Arch(io).outfile = fullfile(pth, fl);
    Arch(io).avg = [];
    Arch(io).ncid = [];
    Arch(io).vid = [];
    
    Arch(io).iens = In.iens;

end
   
    
  
      



