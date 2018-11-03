function varargout = WCVIE2E_physicalmodel(varargin)

%PHYSICALMODEL runs a 2D box model based on Ianson and Allen (2002)

%physicalmodel (outputfile)
%physicalmodel (outputfile, param1, val1, param2, val2...)
%Input = physicalmodel (...)
%Stop=physicalmodel(..., 'stopafterinit',true)

%Input variables:
%
%   REQUIRED:
%   ---------
%
%   outputfile: string, base name of output folder for simulation(s).  The
%               folder name may be modified by the outputextension option
%               (see archiving options, below).  The output folder will
%               contain, at minimum, a dimensions.nc file with temporal and
%               box details and sim0001.nc file with all other
%               output variables.  If ensemble options are used (see
%               archiving options, below, as well as runmixedlayer.m), more
%               simXXXX.nc files may be added, or replaced with
%               post-processed variable files.

%   MODEL GRID:
%   -----------
%
%   dz:         A 2D nz*nx array, the thickness of each box in the vertical
%               dimension
%   
%   dx:         a 2D nz*nx array, the width of each box in the horizontal
%               dimension        
% 
%   dt:         the model time step (days) [1]
%
%   syear:      starting year for simulation (will start on Jan 1 of this
%               year), or 1 x 3 date vector of starting date [1992]
%
%   eyear:      ending year for simulation (will end on Dec 31),
%               or 1 x 3 date vector of ending date [2017]
%
%
%   PHYSICAL PARAMETERS:
%  ---------------------
%   krad1:      the attenuation coefficient (m^-1, value should be
%               positive) for visible radiation (~between 350 nm and 700 nm
%               wavelength).  This is approximately equivalent to the
%               photosynthetically available radiation (PAR). [0.15]
%
%   prad1:      the fraction of incoming solar radiation that falls into
%               visible wavelengths ~ PAR.  By default, it is assumed that
%               roughly 45% of the incoming solar radiation falls into this
%               category (Baker and Frouin, 1987, L&O, 32:6, pp.
%               1370-1377). [0.45]
%
%   krad2:      the attenuation coefficient (m^-1, value should be
%               positive) for non-visible (mainly infra-red) solar
%               radiation.  Water absorbs this radiation very quickly.
%               [1.67]
%
%   alb:        the albedo, or the fraction of incoming radiation reflected
%               from the sea surface.  The default value is 0.079, which is
%               typical for 43 N latitude. (Payne, R.E., 1972, Journal of
%               Atmospheric Sciences, 29:5, pp. 959-969). *If your heat
%               forcing was measured below the water surface, set the
%               albedo to 0. [0.079]
%
%   Lat:        latitude where simulation takes place (used for Coriolis
%               force calculations) [45]
%
%   whgt:       elevation above sea level where wind forcing data was
%               measured (m). [10]  
%
%
%   EXTERNAL FORCING:
%   -----------------
%
%   Qi_input:       Incident solar radiation forcing data. This is an 
%                   (n+1) x (3+nx) matrix. Columns 1-3 hold Year, Month, 
%                   Day and row 1 of the array hold horizontal box names.
%                   Columns 1-3 in row 1 are just placeholders and will be
%                   ignored. The remaining cells hold data for the given 
%                   times and vertical boxes. Radiation is in W.m-2. If all
%                   the forcing data is from the same year, the data will 
%                   be treated as climatological data and will be repeated 
%                   for every simulated year. By default, a dataset 
%                   representing 1976 observations on the Scotia Shelf can 
%                   be used.
%
%
%   airtmp_input:   Air temperature forcing data. Same format as Qi_input.
%                   Temperatures are in deg C. By default, a dataset
%                   representing 1976 observations on the Scotia Shelf can
%                   be used.
%
%   dewptT_input:   Dew point temperature forcing data. Same format as 
%                   Qi_input. Temperatures are in deg C. By default, a 
%                   dataset representing 1976 observations on the Scotia 
%                   Shelf can be used.
%
%   uWndSpd_input:  E-W wind speed forcing data. Same format as Qi_input. 
%                   Speeds in m.d-1. By default, a dataset representing 1976
%                   observations on the Scotia Shelf can be used.
%
%   vWndSpd_input:  N-S wind speed forcing data. Same format as Qi_input. 
%                   Speeds in m.d-1. By default, a dataset representing 1976
%                   observations on the Scotia Shelf can be used.
%
%   wndl_input:     Local longshore wind stress data. This is an nt x 4
%                   matrix, with columns representing year, month, day and
%                   local wind stress (N/m2). Local wind stress is used to
%                   calculate Ekman upwelling.
%
%   wndr_input:     Remote wind forcing data. This is an nt x 4 matrix,
%                   with columns representing year, month, day and remote
%                   wind forcing data (upwelling index) (m.day-1).
%
%   mld_input:      Mixed layer depth (MLD) forcing data (m). Same format 
%                   as Qi_input. 
%
%   entrnmnt_input: entrainment rate data (day-1). Same format 
%                   as Qi_input. 
%
%   INITIAL CONDITIONS:
%   -------------------
%        
%   t_input:        Initial temperature profiles for simulation.
%                   Data is an nz x nx matrix with initial temperatures.
%              
%   s_input:        Initial salinity profiles for simulation.
%                   Data is an nz x nx matrix with initial salinities.
%
%   PHYSICAL LATERAL BOUNDARY DATA:
%   -------------------------------
%   LBtmp_input:    Lateral boundary temperatures. Columns 1-3 of this array
%                   hold a year, month and day corresponding to the dates
%                   of the boundary conditions, and row 1 of the array
%                   holds the origin: Column 1 = open ocean upper
%                   layer; Column 2 = open ocean lower layer; Column 3 = 
%                   rain; Column 4 = freshwater from run-offs; 
%                   Column 5 = VICC. Columns 1-3 in row 1 are 
%                   just placeholders and will be ignored. The remaining cells
%                   hold lateral boundary temperature data (deg C) for the
%                   given times and sources. If all the boundary data
%                   is from the same year, the data will be treated as
%                   climatological data and will be repeated for every
%                   simulated year.
%
%   LBsal_input:    Lateral boundary salinities. Format is the same as for
%                   LBtmp_input, with salinity data in ppt. If all the
%                   relaxation data is from the same year, the data will be
%                   treated as climatological data and will be repeated for
%                   every simulated year.
%
%   PHYSICAL RELAXATION DATA:
%   -------------------------
%
%   srelax:         Relaxtion data for salinity. This is an (nz x nx x n) 
%                   x 6 cell array. Columns 1-6 of this cell array hold a 
%                   year, a month, a day, the box in the x dimension, the 
%                   box in the z dimension and salinity data(ppt)for the 
%                   given time and boxes towards which the modeled salinity
%                   will be relaxed. If all the relaxation data is from the
%                   same year, the data will be treated as climatological 
%                   data and will be repeated for every simulated year. 
%                   If empty, no relaxation will be done.
%
%   srelaxtime:     Timescale for salinity relaxation (days) [i.e. 30 d]
%                   
%
%   trelax:         Relaxation data for temperature.  Format is the same as 
%                   for srelax, with temperature data in deg C. If all the
%                   relaxation data is from the same year, the data will be
%                   treated as climatological data and will be repeated for
%                   every simulated year.  If empty, no relaxation will be
%                   done. 
%
%   trelaxtime:     Timescale for temperature relaxation (days) [i.e. 30 days]
%                   
%
%   ARCHIVING:
%   ----------
%
%   tarch:          the archiving interval (days). Data is averaged over
%                   tarch days.
%                   If tarch isn't scalar, then multiple output folders will be
%                   created; beginarchive and endarchive must be the same size
%                   as tarch. [30]
%
%   beginarchive:   datenumber indicating what model date to begin recording
%                   to an output file, or NaN to indicate that archiving begins
%                   immediately. [NaN]
%
%   endarchive:     datenumber indicating what model date to stop recording
%                   output to an output file, or NaN to indicate archiving until
%                   the simulation ends. [NaN]
%
%   outputextension:cell array of strings, same size as tarch, beginarchive, 
%                   and endarchive.  This string is appended to the 
%                   outputfile string if multiple output folders are 
%                   indicated by the other archiving variables.  If empty 
%                   and tarch is non-scalar,then numeric suffixes will be 
%                   used (i.e. outputfile1, outputfile2, etc.)
%
%   nens:           number of ensemble members in a set of mixed_layer runs, to
%                   be saved to the same output folder.  All ensemble member
%                   runs *must* use the same spatial and temporal grid,
%                   archiving options, and biological module; all other
%                   parameters can be varied between simulations. [1]
%
%   iens:           index of the current ensemble member.  Determines the
%                   number assigned to the output file name, and concatenation
%                   order if postprocessing is used (see runmixedlayer.m) [1]
%
%   stopafterinit: logical scalar.  If true, simulation is terminated after
%                  the initialization process, and no forward integration is
%                  performed.  This is for debugging purposes.  Using this
%                  option also leads to different output if mixed_layer is
%                  called with an output variable (see below) [false]
%   BIOLOGY:
%   --------
%
%   biofun:         Function handle to biological module.  If empty, the model
%                   will run without any biology. [empty] 
%
%   *var*relax:     Relaxation data for any biological state variable, where
%                   *var* corresponds to the short name of the variable.
%                   Format is the same as for srelax.  If not included or
%                   empty, no relaxation will be done.
%
%   *var*flux:      Additional flux into (or out of) any biological state
%                   variable, where *var* corresponds to the short name of the
%                   variable.  Format is the same as for srelax, in units of
%                   stateVariableUnit/s.
%
%   brelaxtime:     Timescale for relaxation of biological state variables. (s)
%                   [2592000]
%
%   openbottom:     Logical scalar. The value changes the way the biological
%                   tracer variables interact in the bottom cell.  If false, a
%                   no-flux condition is set at the bottom; mixing and vertical
%                   movement will be conservative for biological tracers.  If
%                   true, the values of biological variables are held constant
%                   in the bottom cell, and material sinks through the
%                   bottom boundary; mass is not conserved under these
%                   conditions. [false]
%
%   BIOLOGICAL RELAXATION DATA:
%   ---------------------------
%
%
%   OTHER
%   -----
%
%   verbose:        logical scalar.  If true, progress statements will be
%                   printed to the screen.  If false, nothing will be printed.
%                   [true]
%
%   RESTART
%   -------
%
%   Note: these options are a bit kludgy; hot starts should only be run
%   when just changing the forcing, but keeping most input (specifically,
%   archiving variables) the same.
%
%   hotstartdn:     scalar datenumber, indicating time point at which
%                   conditions should be saved, for use in restarting a
%                   simulation at that point []
%
%   hotstartsave:   name of .mat file where hot start data will be saved
%                   (.mat extension should be included).  Only valid when a
%                   value is supplied for hotstartdn. ['mlhotstart.mat']
%
%   hotstartload:   name of file to use when restarting.  Simulation will
%                   start at the last time index in the file.  All output will
%                   be the same as it would have been had the simulation been
%                   run from the beginning.
%
% Output variables:
%
%   Input:          1 x 1 structure holding the input variables used for the
%                   run.  Returning this variable allows you to see all inputs
%                   used, including those set by defaults.
%
%   Stop:           1 x 1 structure holding the variables used within a
%                   mixed_layer run.  This is returned only if the
%                   'stopafterinit' flag is set to true, and is useful for
%                   debugging purposes.

% Copyright 2011-2015 Kelly Kearney, Charlie Stock
% kakearney@gmail.com, charles.stock@noaa.gov

%--------------------------
% Setup
%--------------------------
%
% Check input

In = WCVIE2E_parseinput(varargin{:});

if isfield(In, 'stopafterinit') % flag to stop after initialization
    stopafterinit = true;
    In = rmfield(In, 'stopafterinit');
else
    stopafterinit = false;
end

if In.verbose
    fprintf('\n--------------------\n2D Physical Model\n--------------------\n\n');
    fprintf('%s: %d\n\n', In.outputfile, In.iens);
end

%--------------------------
% Initialization
%--------------------------

if In.verbose
    fprintf('Initializing model...\n');
end

% Set constant parameters and initial conditions for simulation

[Grd, Ht, Wnd, VMxng, TS, LB, Arch] = WCVIE2E_initialize(In);

% If biology is included, set initial conditions and extra parameters for
% biological calculations

% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%

if In.hasbio;
    [Bio.bio, Bio.ismixed, Bio.bottomval, Bio.Params, Bio.names, ...
        Bio.diagnames] = feval(In.biofun, 'init', In, Grd);
    Bio.ndiag = size(Bio.diagnames,1);
    if Bio.ndiag > 0
        Bio.diag = zeros(Grd.nz,Bio.ndiag); % placeholder KAK note: Charlie's
        %new version runs biofun here to get ititial conditions
    end
    
%     Con = bioconservecheck(In.outputfile);
%     Con.dz = -diff(Grd.zp);
end

% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%

% Set up relaxation of biological variables (and additional outside fluxes)

if In.hasbio
    Bio = initbiorelax(Grd, Bio, In);
end
% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%

% Stop-after-init (for troubleshooting purposes)

if stopafterinit
    stoptempfile = tempname;
    save(stoptempfile);         % Easiest way to gather all variables into 
    A = load(stoptempfile);     % a structure is to save and reload
    delete([stoptempfile '.mat']);
    varargout{1} = A;
    if In.verbose
       fprintf('Initialization-only run; exiting mixed_layer\n');
    end
    return
end

% Set up print strings for screen display

simdatestr = datestr(datenum(Grd.start_date) + Grd.time, '  mm/dd/yyyy\n');
dummystr = [repmat(' ', 1, size(simdatestr,2))-1 '\n'];
erasestr = repmat('\b', 1, size(simdatestr,2));

%--------------------------
% Mixing simulation
%--------------------------

mltimer = tic;

if ~isempty(In.hotstartload)
    Tmp = load(In.hotstartload);
    
    for ia = 1:length(Arch)
        Arch(ia).avg = Tmp.Arch(ia).avg;
    end
    
    TS.T = Tmp.TS.T;
    TS.S = Tmp.TS.S;
    TS.Sig = Tmp.TS.Sig;
    
 % --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%
 
    if In.hasbio
        Bio.bio = Tmp.Bio.bio;
        Bio.diag = Tmp.Bio.diag;
    end
% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%
    
    tidx = (Tmp.Grd.savehot+1):Grd.nt;
else
    tidx = 1:Grd.nt;
end    

for it=tidx
    
%----------------------------------------------------
% Extract forcing data, boundary data and relaxation
% data at current time step.
%----------------------------------------------------

Ht.qi=Ht.Qi.data(it,:);% row vector size 2
Ht.airtmp=Ht.Airtmp.data(it,:);% row vector size 2
Ht.dewptT=Ht.DewptT.data(it,:);% row vector size 2
Ht.qo=Ht.Qo.data(it);%scalar
Ht.meanqi=Ht.MeanQi.data(it,:);% row vector size 2. not sure about the dim of this one as 'smooth' was used

Wnd.uwndspd=Wnd.UWndSpd.data(it,:);% row vector size 2
Wnd.vwndspd=Wnd.VWndSpd.data(it,:);% row vector size 2
Wnd.wndspd10=Wnd.WndSpd10.data(it,:);% row vector size 2
Wnd.wndl=Wnd.Wndl.data(it);%scalar
Wnd.wndr=Wnd.Wndr.data(it);%scalar

VMxng.mld=VMxng.Mld.data(it,:);% row vector size 2
VMxng.entrnmnt=VMxng.Entrnmnt.data(it,:);% row vector size 2

LB.physics.lbtmp=LB.physics.LBtmp.data(it,:);%row vector size 5
LB.physics.lbsal=LB.physics.LBsal.data(it,:);%row vector size 5

if In.hassrelax
    TS.srelax= TS.Srelax.data(:,:,it);% 3 x 2 matrix
end

if In.hastrelax
    TS.trelax= TS.Trelax.data(:,:,it);% 3 x 2 matrix
end

% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%
if In.hasbio 
    nb = size(Bio.names,1);
        Bio.brelax = cell(nb,1);
        Bio.extraflx = cell(nb,1);
        for ib = 1:nb
            if Bio.hasrelax(ib)
                Bio.brelax{ib} = Bio.Relax(ib).data(:,:,it);
            end
            if Bio.hasflux(ib)
            Bio.extraflx{ib} = Bio.ExtraFlux(ib).data(:,:,it);
            end
        end
end
% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%
 
%--------------------------
% Calculate some heat and
% wind parameters
%--------------------------

% Calculate incident, sensible, latent, and longwave heat fluxes, based
% on heat forcing/stratification at the start of each time step
     
[Ht.qi, Ht.qs, Ht.ql, Ht.qlw] = WCVIE2E_calcheat(Ht.qi, ...
Ht.airtmp, Ht.dewptT, TS.T(1,:), Wnd.wndspd10, ...
Ht.qo, Ht.meanqi, In.alb, Grd.nx);

% Incident heat flux is distributed throughout the water column. prad
% is the fraction that is photosynthetically-active, attenuated
% according to krad1.  The remaining fraction is attenuated according
% to krad2.  Qi is also adjusted downward based on albedo.

for i = 1:Grd.nx
Ht.solhflx(:,i) = -(1-In.alb).* Ht.qi(i).*diff(...
        In.prad1.*exp(In.krad1.*Grd.zp(:,i)) + (1-In.prad1).*exp(In.krad2.*Grd.zp(:,i)));
end

% Sensible, latent, and longwave fluxes are applied at the surface
    
Ht.srfhflx = Ht.qs + Ht.ql + Ht.qlw;

% Heat fluxes (W/m^2) are converted to temperature fluxes (deg C/s)
    
Cp = 4125; % Typical specific heat for water (J kg^-1 degC^-1) @10 deg., 32 ppt)
    
Ht.sol_Tflx = Ht.solhflx./(TS.Sig.*Cp.*In.dz);
Ht.srf_Tflx = Ht.srfhflx./(TS.Sig(1,:).*Cp.*In.dz(1,:));

%--------------------------
% Archiving process
%--------------------------
                      
% Archived variables description: variable, short name, long name,
% units

  datatoarchive = {...
        TS.T                        'temp'      'temperature'                               'deg C'
        TS.S                        'sal'       'salinity'                                  'psu'
        TS.Sig                      'sig'       'density'                                   'kg m^-3'
        Ht.qo                       'Qo'        'clear sky irradiance'                      'W m^-2'
        Ht.meanqi                   'mQi'       'mean observed daily irradiance'            'W m^-2'
        Ht.qi                       'Qi'        'incident heat flux'                        'W m^-2'
        Ht.qs                       'Qs'        'sensible heat flux'                        'W m^-2'
        Ht.ql                       'Ql'        'latent heat flux'                          'W m^-2'
        Ht.qlw                      'Qlw'       'longwave heat flux'                        'W m^-2'
        Ht.solhflx                  'Qsol'      'incident heat flux (water column)'         'W m^-2'
        Ht.srfhflx                  'Qsrf'      'surface heat flux'                         'W m^-2'
        Wnd.wndspd10                'wndspd10'  'speed velocity 10 m above sea level'       'm d^-1'
        };       
    
    % If running biology, add state variables and diagnostics to this
    % matrix
    
% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%

    if In.hasbio
        if Bio.ndiag > 0
            biodata = num2cell([Bio.bio Bio.diag],1)';
            datatoarchive = [datatoarchive; [biodata [Bio.names; Bio.diagnames]]];
        else
            biodata = mat2cell(Bio.bio, size(Bio.bio,1), ones(1, size(Bio.bio,2)))';
            datatoarchive = [datatoarchive; [biodata Bio.names]];
        end
    end    
% --------------TO VERIFY WHEN WORKING ON THE BIOLOGICAL MODULES---------%
    
     % Archiving
    
    for io = 1:length(Arch)
        [Arch(io).avg, Arch(io).ncid, Arch(io).vid] = WCVIE2E_archivemldata(Grd, Arch(io), it, datatoarchive, tidx(1));
    end
    if it == tidx(1)
        cu = onCleanup(@() closefiles([Arch.ncid]));
    end
    
    % Print simulation date to screen
    
    if it == tidx(1)
		if In.verbose
	        fprintf('Running simulation...\n');
	        fprintf(dummystr);
		end 
    end
    
    if Grd.newday(it) 
		if In.verbose
	        fprintf([erasestr simdatestr(it,:)]);
		end
    end  
    
%--------------------------
% Calculate mixing, 
% sinking, and growth over 
% this time step
%--------------------------
    
% Mix temperature and salinity     
    
 TS.T = mixtracer(TS.T, In, Grd, VMxng.mld, VMxng.entrnmnt, ...
     LB.physics.lbtmp, Wnd.wndl, Wnd.Wndr, TS.Sig, it,   ...
                     'sflux', Ht.srf_Tflx, 'bflux', 0, ...
                     'source', Ht.sol_Tflx);
                 
 TS.S = mixtracer(TS.S, In, Grd, VMxng.mld, VMxng.entrnmnt, ...
     LB.physics.lbsal, Wnd.wndl, Wnd.Wndr, TS.Sig, it,   ...
                     'sflux', 0, 'bflux', 0); 
    
    


% Relax temperature and salinity

    if In.hastrelax
        TS.T = TS.T + (TS.trelax - TS.T)*In.dt./In.trelaxtime;%trelaxtime=1?
    end
    
    if In.hassrelax
        TS.S = TS.S + (TS.srelax - TS.S)*In.dt./In.srelaxtime;%trelaxtime=1?
    end
    
 %----------------------------------------
 %TO CHECK WHEN GOING THROUGH THE BIOLOGY
 %----------------------------------------  
 % Relax biology
    
    if In.hasbio
        for ib = 1:nb
            if Bio.hasrelax(ib)
                leavealone = isnan(Bio.brelax{ib});
                Bio.bio(~leavealone,ib) = Bio.bio(~leavealone,ib) + (Bio.brelax{ib}(~leavealone) - Bio.bio(~leavealone,ib))*In.dt/In.brelaxtime;
            end
        end
    end
 %----------------------------------------
 %TO CHECK WHEN GOING THROUGH THE BIOLOGY
 %---------------------------------------- 
 
 % Calculate density

  TS.Sig = WCVIE2E_sw_dens0(TS.S, TS.T);
  
 %----------------------------------------
 %TO CHECK WHEN GOING THROUGH THE BIOLOGY
 %---------------------------------------- 
  
  % Mix biological variables

    if In.hasbio
        
        % Turbulent mixing
        
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 1);
            
        premixbio = Bio.bio;
        for ibio = 1:size(Bio.bio,2)
            if Bio.ismixed(ibio)
                if isnan(Bio.bottomval(ibio))
                    Bio.bio(:,ibio) = mixtracer(Bio.bio(:,ibio), ...
                                    Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
                                    'bflux', 0, 'sflux', 0);
                else
                    Bio.bio(:,ibio) = mixtracer(Bio.bio(:,ibio), ...
                                    Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
                                    'bval', Bio.bottomval(ibio), ...
                                    'sflux', 0);
                end
            end
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 2);
            
        % NaN check
        
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
        % Open bottom indicator (I do this here rather than during initial
        % setup due to some issues that arose in my bio models when state
        % variables were 0 at the bottom.  By waiting until after one
        % mixing calculation, I ensure that all bio has at least a
        % minuscule amount of biomass.  A bit of a hack, and may not work
        % in all cases, so might need to look into this further later).

        if it == 1 && In.openbottom
            Bio.bottomval = Bio.bio(end,:);
        end
        
        % Other vertical movement
        
        PhysParams = Ts;
        PhysParams.par = Ht.Qi.*In.prad1;
        PhysParams.par24 = Ht.meanQi(it).*In.prad1;
        PhysParams.kpar = In.krad1;
        
        GrdParams = struct('z', Grd.z, 'dz', In.dz, 't', Grd.time(it), 'dt', In.dt);
                
        Bio.wsink = feval(In.biofun, 'vertmove', Bio.bio, ...
                              PhysParams, Bio.Params, GrdParams);
            
        presinkbio = Bio.bio;
        for ibio = 1:size(Bio.bio,2)
            Bio.bio(:,ibio) = verticalflux(Bio.bio(:,ibio), Bio.wsink(:,ibio), In.dt, In.dz, In.openbottom);
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 3);
     
          
    end

 %----------------------------------------
 %TO CHECK WHEN GOING THROUGH THE BIOLOGY
 %---------------------------------------- 

% Solve for biological sources/sinks and step forward

    if In.hasbio
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
%         PhysParams = Ts;
%         PhysParams.par = Ht.Qi.*In.prad1;
%         PhysParams.par24 = Ht.meanQi(it).*In.prad1;
%         PhysParams.krad1 = In.krad1;
%         
%         GrdParams = struct('z', Grd.z, 'dz', In.dz, 't', Grd.time(it), 'dt', In.dt);
%         
        [Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
                              PhysParams, Bio.Params, GrdParams);
        
%         [Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
%                                   Ht.meanQi(it).*In.prad1, Ts.T, Grd.z, In.dz, ...
%                                   Bio.Params, Grd.time(it), In.dt);
                              
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 4);
    end
    
   if In.hasbio
        for ib = 1:nb
            if Bio.hasflux(ib)
                dbdtExtra = interp2(Bio.ExtraFlux(ib).t, Bio.ExtraFlux(ib).z, Bio.ExtraFlux(ib).data, Grd.time(it), Grd.z);
                Bio.bio(:,ib) = Bio.bio(:,ib) + dbdtExtra .* In.dt;
            end
        end
    end

   if it == Grd.savehot
        save(In.hotstartsave);
    end

end
 
runtime = toc(mltimer);

% Print loop runtime

if In.verbose
	fprintf('\nMixed layer model completed successfully: %f s\n', runtime);
end

% Return input variables if output requested

if nargout == 1
    varargout{1} = In;
end

%--------------------------
% Cleanup
%--------------------------

% Debugging

if In.hasbio && isfield(Bio.Params, 'cfid')
    fclose(Bio.Params.cfid);
end


function closefiles(id)
for ii = 1:length(id)
    netcdf.close(id(ii));
end

