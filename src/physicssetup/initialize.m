function [Grd, TS, Arch]=initialize(In)
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
%   Arch:   structure holding variables related to archive (i.e. output).
%           The archiving period refers to the In.tarch input.
%
%           startdate:  nbin x 6 array or nbin x 1 array, date vectors or 
%                       date numbers corresponding to
%                       start of each archiving time step.
%
%           enddate:    nbin x 6 array or nbin x 1 array, date vectors or
%                       date numbers corresponding to
%                       end of each archiving time step.
%
%           middate:    nbin x 6 array or nbin x 1 array, date vectors or
%                       date numbers corresponding to
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

Grd.tmax = (dnend - dnstart)*86400; %number of seconds elapsed from 1 Jan 1992 at midnight to 1 jan 2017 at midnight
Grd.nt = floor(Grd.tmax/In.dt); %total number of time steps throughout the entire simulation

Grd.time=(0:Grd.nt)*In.dt;

% For screen print counter: print progress at the end of each day

daydiff=diff(floor(Grd.time/86400+datenum(Grd.start_date)));
Grd.newday=[true logical(daydiff)];

% time step properties for interpolating forcing datasets for ODE solver

Grd.datant=floor(Grd.tmax./In.datadt);
Grd.datatime=(0:Grd.datant).*In.datadt;

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

TS.Sig = sw_dens0(TS.S, TS.T);

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
    
    %binedge{io} = 0:In.tarch(io):max(Grd.time);    
    dn = dnstart + Grd.time./86400;
    datevector = datevec(dn);
    isedge = ismember(datevector(:,3:6),[1 0 0 0],'rows');
    binedge{io} = Grd.time(isedge);    
    
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
    
    % Dates
    
    Arch(io).dateedge = dnstart + binedge{io}./86400;
    %Arch(io).dateedge = dnstart + Grd.time./86400;
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
    Arch(io).file = [];
    Arch(io).ncid = [];
    Arch(io).vid = [];
    
    Arch(io).iens = In.iens;

end
   

    

      



