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
%   dt:         the model time step (seconds) [21600 ~ 6 hrs]
%
%   datadt:     Time-step required for interpolating forcing data. A finer
%               time-step might be required when using certain ODE solvers.
%               [1800 ~ 0.5 hrs]
%               
%   syear:      starting year for simulation (will start on Jan 1 of this
%               year), or 1 x 6 date vector of starting date [1992]
%
%   eyear:      ending year for simulation (will end on Dec 31),
%               or 1 x 6 date vector of ending date [2017]
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
%               force calculations) [48]
%
%   whgt:       elevation above sea level where wind forcing data was
%               measured (m). [10]  
%
%
%   PHYSICAL EXTERNAL FORCING:
%   -------------------------
%
%   Qi_input:       Incident solar radiation forcing data (i.e. downward 
%                   shortwave radiation fluxes). This is an (n+1) x (6+nx) 
%                   matrix. Columns 1-6 hold Year, Month, Day, Hour, Minute
%                   and Seconds and the rest of the columns(7:end) holds 
%                   data for each horizontal box. Row 1 of the array hold 
%                   date time vectors (Y, M, D, H, MIN, S) and horizontal 
%                   box names. Row 1 will be ignored. The remaining cells 
%                   hold data for the given times and horizontal boxes. 
%                   Radiation is in W.m-2. If all the forcing data is from 
%                   the same year, the data will be treated as climatological 
%                   data and will be repeated for every simulated year.
%
%   airtmp_input:   Air temperature forcing data. Same format as Qi_input.
%                   Temperatures are in degrees. If all the forcing data is from 
%                   the same year, the data will be treated as climatological 
%                   data and will be repeated for every simulated year.
%
%   dewptT_input:   Dew point temperature forcing data. Same format as
%                   Qi_input. Temperatures are in degrees. If all the forcing
%                   data is from the same year, the data will be treated as 
%                   climatological data and will be repeated for every 
%                   simulated year.
%
%   uWndSpd_input:  E-W wind speed forcing data. Same format as Qi_input. 
%                   Speeds in m.s-1. 
%
%   vWndSpd_input:  N-S wind speed forcing data. Same format as Qi_input. 
%                   Speeds in m.s-1. 
%
%   mld_input:      Mixed layer depth (MLD) forcing data (m). Same format 
%                   as Qi_input. Remove this dataset once I get the code
%                   from Ze.
%
%  entrnmnt_input:  entrainment rate data (second-1). Same format as Qi_input. 
%                   Remove this dataset once I get the code from Ze.
%
%  p_input       :  precipitation rate data (m.s^-1). Same format as
%                   Qi_input.
%
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
%   LBtmp_input:    Lateral boundary temperatures. Columns 1-6 of this array
%                   hold a year, month, day, hour, minute and second 
%                   corresponding to the dates of the boundary conditions, 
%                   and the remaining columns hold lateral boundary 
%                   temperature data (deg C) for the given times and sources. 
%                   Column 7 = open ocean upper layer; Column 8 = open ocean lower layer; 
%                   Column 9 = rain over the shelf; Column 10= rain over the slope; 
%                   Column 11 = freshwater from run-offs;Column 12 = VICC. 
%                   Row 1 will be ignored. If all the boundary data is from 
%                   the same year, the data will be treated as climatological 
%                   data and will be repeated for every simulated year.
%
%   LBsal_input:    Lateral boundary salinities. Format is the same as for
%                   LBtmp_input, with salinity data in ppt. If all the
%                   relaxation data is from the same year, the data will be
%                   treated as climatological data and will be repeated for
%                   every simulated year.
%
%   DIAGNOSTIC:
%   -----------
%
%   Tdiag:         logical scalar specifying whether a diagnostic should be
%                  run for temperature. It basically writes out output
%                  variables that directly affect temperature for every time
%                  step. [false]
%
%   Sdiag:         logical scalar specifying whether a diagnostic should be
%                  run for salinity. It basically writes out output
%                  variables that directly affect salinity for every time
%                  step. [false]
%
%   ARCHIVING:
%   ----------
%
%   tarch:          the archiving interval (seconds). Data is averaged over
%                   tarch seconds.
%                   If tarch isn't scalar, then multiple output folders will be
%                   created; beginarchive and endarchive must be the same size
%                   as tarch. [2592000]
%
%   beginarchive:   datenumber indicating what model date to begin recording
%                   to an output file, or NaN to indicate that archiving begins
%                   immediately. [NaN]
%
%   endarchive:     datenumber indicating what model date to stop recording
%                   output to an output file, or NaN to indicate archiving until
%                   the simulation ends. [NaN]
%
%   spatialarch:    logical scalar indicating whether outputs are
%                   spatially-resolved [true]
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
%   openbottom:     Logical scalar. The value changes the way the biological
%                   tracer variables interact in the bottom cell.  If false, a
%                   no-flux condition is set at the bottom; mixing and vertical
%                   movement will be conservative for biological tracers.  If
%                   true, the values of biological variables are held constant
%                   in the bottom cell, and material sinks through the
%                   bottom boundary; mass is not conserved under these
%                   conditions. [false]
%
%   OTHER
%   -----
%
%   verbose:        logical scalar.  If true, progress statements will be
%                   printed to the screen.  If false, nothing will be printed.
%                   [true]
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

mltimer = tic;

%--------------------------------------------------------------------------
% Setup
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------

if In.verbose
    fprintf('Initializing model...\n');
end

% Set forcing datasets, lateral boundary conditions datasets and initial 
% conditions for simulation

[Grd, Ht, Wnd, VMxng, TS, LB, P, Arch, Tdiag, Sdiag] = WCVIE2E_initialize(In);

% If biology is included, set initial conditions and extra parameters for
% biological calculations

if In.hasbio
    [Bio.bioin, Bio.bio, Bio.ismixed, Bio.bottomval, Bio.Params, Bio.names, ...
        Bio.diagnames, Bio.bgcdata, Bio.biots, Bio.LB] = feval(In.biofun, 'init', In, Grd);
    Bio.ndiag = size(Bio.diagnames,1);
    if Bio.ndiag > 0
        Bio.diag = zeros(Grd.nz,Grd.nx,Bio.ndiag); % placeholder KAK note: Charlie's new version runs biofun here to get ititial conditions
    end
    
%     Con = bioconservecheck(In.outputfile);
%     Con.dz = -diff(Grd.zp);
end

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

simdatestr = datestr(datenum(Grd.start_date) + Grd.time./86400, '  mm/dd/yyyy\n');
dummystr = [repmat(' ', 1, size(simdatestr,2))-1 '\n'];
erasestr = repmat('\b', 1, size(simdatestr,2));


tidx = 1:Grd.nt;

for it=tidx
    
%--------------------------------------------------------------------------
% Archiving process
%--------------------------------------------------------------------------

% Archived variables description: variable, short name, long name,
% units

%datatoarchive = {...
% TS.T                        'temp'      'temperature'                               'deg C'
% TS.S                        'sal'       'salinity'                                  'psu'
% TS.Sig                      'sig'       'density'                                   'kg m^-3'
% };    

% If running biology, add state variables and diagnostics to this
% matrix
    
if In.hasbio
   if Bio.ndiag > 0
      biodata = num2cell(cat(3,Bio.bio,Bio.diag),[1 2]);
      biodata = reshape(biodata,size(biodata,3),[],[]);
      datatoarchive = [datatoarchive; [biodata [Bio.names; Bio.diagnames]]];
   else
      biodata = mat2cell(Bio.bio, size(Bio.bio,1), size(Bio.bio,2),...
          ones(1, size(Bio.bio,3)));
      biodata = reshape(biodata,size(biodata,3),[],[]);
      datatoarchive = [datatoarchive; [biodata Bio.names]];
   end
end
               
% Archiving
    
%for io = 1:length(Arch)
%    [Arch(io).avg, Arch(io).file, Arch(io).ncid, Arch(io).vid] = WCVIE2E_archivemldata(Grd, Arch(io), it, datatoarchive, tidx(1));
%end

%if it == tidx(1)
%   cu = onCleanup(@() closefiles([Arch.ncid]));
%end


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

%--------------------------------------------------------------------------
% Stepping forward- Resolving ODE for temperature and salinity
%--------------------------------------------------------------------------

%calculate new temperature and salinity

[TS.T, rk1_T, V1T, H1T, X1T, rk2_T, V2T, H2T, X2T, rk3_T, V3T, H3T, X3T,...
rk4_T, V4T, H4T, X4T, srf_Tflx1, sol_Tflx1, srf_Tflx2, sol_Tflx2, ...
srf_Tflx3, sol_Tflx3, srf_Tflx4, sol_Tflx4] = rk4(@mixtracerTS_ode, it, ...
TS.T, In, Grd, VMxng.Mld, VMxng.Entrnmnt, LB.tmp, TS.Sig, Wnd.Tauy2, P, ...
'sflux', Ht.Qi, Ht.airtmp, Ht.dewptT, Wnd.Wspd10, Ht.Qo,'source',Ht.Qi);
            
[TS.S, rk1_S, V1S, H1S, X1S, rk2_S, V2S, H2S, X2S, rk3_S, V3S, H3S, X3S,...
    rk4_S, V4S, H4S, X4S] = rk4(@mixtracerTS_ode, it, TS.S, In, Grd, ...
    VMxng.Mld, VMxng.Entrnmnt,LB.sal, TS.Sig, Wnd.Tauy2, P);

% Mix biological variables

if In.hasbio
        
% Turbulent mixing
        
if any(isnan(Bio.bio(:)))
   warning('ML: NaN in biology')
end
        
%Con = bioconservecheck(Con, Bio.bio, it, 1);
            
 for ibio = 1:size(Bio.bio,3)
     if Bio.ismixed(ibio)
        if isnan(Bio.bottomval{ibio}.data)%find a solution to add function all
           Bio.bio(:,:,ibio) = rk4(@mixtracerBio_ode, it, Bio.bio(:,:,ibio), ...
                               In, Grd, VMxng.Mld, VMxng.Entrnmnt,...
                               Bio.LB{ibio}, TS.Sig, Wnd.Tauy, P); % change lateral boundary conditions with bio  
        else
           Bio.bio(:,:,ibio) = rk4(@mixtracerBio_ode, it, Bio.bio(:,:,ibio), ...
                               In, Grd, VMxng.Mld, VMxng.Entrnmnt,...
                               Bio.LB{ibio}, TS.Sig, Wnd.Tauy, P,...
                               'bval', Bio.bottomval{ibio});% change lateral boundary conditions with bio  
        end
     end
 end
 
end
        
% Calculate new density

TS.Sig = WCVIE2E_sw_dens0(TS.S, TS.T);


if In.hasbio

%Con = bioconservecheck(Con, Bio.bio, it, 2);
            
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
    for i=1:size(Bio.bio,3)
   Bio.bottomval{i}.data = repmat(Bio.bio(end,:,i),Grd.datant+1,1);
    end
end
    
% set parameter structures for biology

PhysParams = TS;
PhysParams.par24=Ht.Qi.data(Ht.Qi.t==Grd.time(it),:).*In.prad1;%1 x nx array
PhysParams.kpar = In.krad1;

        
GrdParams = struct('z', Grd.z, 'dz', In.dz, 'dx', In.dx, 'area',...
    Bio.bioin.area, 't', Grd.time(it), 'dt', In.dt, 'sdate', Grd.start_date);

% Other vertical movement
                
Bio.wsink = feval(In.biofun, 'vertmove', Bio.Params, GrdParams);
            
for ibio = 1:size(Bio.bio,3)
    Bio.bio(:,:,ibio) = verticalflux(Bio.bio(:,:,ibio), Bio.wsink(:,:,ibio), In.dt, In.dz, In.openbottom);
end
        
               
end

% Solve for biological sources/sinks and step forward

if In.hasbio
   if any(isnan(Bio.bio(:)))
      warning('ML: NaN in biology')
   end
        
%PhysParams = Ts;
%PhysParams.par = Ht.Qi.*In.prad1;
%PhysParams.par24 = Ht.meanQi(it).*In.prad1;
%PhysParams.krad1 = In.krad1;
%         
% GrdParams = struct('z', Grd.z, 'dz', In.dz, 't', Grd.time(it), 'dt', In.dt);
%         
[Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
                            PhysParams, Bio.Params, Bio.bgcdata, ...
                            Bio.biots, GrdParams);
        
%[Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
%Ht.meanQi(it).*In.prad1, Ts.T, Grd.z, In.dz, ...
%Bio.Params, Grd.time(it), In.dt);
                              
if any(isnan(Bio.bio(:)))
   warning('ML: NaN in biology')
end
        
%Con = bioconservecheck(Con, Bio.bio, it, 4);

end

%--------------------------------------------------------------------------
%OPTIONAL : Run some diagnostic tests
%--------------------------------------------------------------------------

% Diagnostic tests for temperature

if In.hasTdiag
     
%Compile all data of interest for diagnostic purposes
    datatodiag = {...
       TS.T                        'temp'      'temperature'                               'deg C'
       TS.Sig                      'sig'       'density'                                   'kg.m-3'
       rk1_T                       'ode1T'     '1st order ode for T'                       'deg C.s^-1'
       V1T                         'odeV1T'    'vertical flux 1st order ode T'             'deg C.s^-1'
       H1T                         'odeH1T'    'horizontal flux 1st order ode T'           'deg C.s^-1'
       X1T                         'odeX1T'    'advection flux 1st order ode T'            'deg C.s^-1'
       srf_Tflx1
       sol_Tflx1
       rk2_T                       'ode2T'     '2nd order ode for T'                       'deg C.s^-1'       
       V2T                         'odeV2T'    'vertical flux 2nd order ode T'             'deg C.s^-1'
       H2T                         'odeH2T'    'horizontal flux 2nd order ode T'           'deg C.s^-1'
       X2T                         'odeX2T'    'advection flux 2nd order ode T'            'deg C.s^-1'       
       srf_Tflx2
       sol_Tflx2
       rk3_T                       'ode3T'     '3rd order ode for T'                       'deg C.s^-1'
       V3T                         'odeV3T'    'vertical flux 3rd order ode T'             'deg C.s^-1'
       H3T                         'odeH3T'    'horizontal flux 3rd order ode T'           'deg C.s^-1'
       X3T                         'odeX3T'    'advection flux 3rd order ode T'            'deg C.s^-1'            
       srf_Tflx3
       sol_Tflx3
       rk4_T                       'ode4T'     '4th order ode for T'                       'deg C.s^-1'
       V4T                         'odeV4T'    'vertical flux 4th order ode T'             'deg C.s^-1'
       H4T                         'odeH4T'    'horizontal flux 4th order ode T'           'deg C.s^-1'
       X4T                         'odeX4T'    'advection flux 4th order ode T'            'deg C.s^-1'             
       srf_Tflx4
       sol_Tflx4
       };

% Write data to file
    
[Tdiag.file, Tdiag.ncid, Tdiag.vid] = diagoutput(Grd, Tdiag, it, datatodiag, tidx(1));

if it == tidx(1)
   cu = onCleanup(@() closefiles([Tdiag.ncid]));
end
   
end

% Diagnostic tests for salinity

if In.hasSdiag
     
%Compile all data of interest for diagnostic purposes
    datatodiag = {...
       TS.S                        'sal'       'salinity'                                  'psu'
       TS.Sig                      'sig'       'density'                                   'kg.m-3'
       rk1_S                       'ode1S'     '1st order ode for S'                       'psu.s^-1'
       V1S                         'odeV1S'    'vertical flux 1st order ode S'             'psu.s^-1'
       H1S                         'odeH1S'    'horizontal flux 1st order ode S'           'psu.s^-1'
       X1S                         'odeX1S'    'advection flux 1st order ode S'            'psu.s^-1'
       rk2_S                       'ode2S'     '2nd order ode for S'                       'psu.s^-1'       
       V2S                         'odeV2S'    'vertical flux 2nd order ode S'             'psu.s^-1'
       H2S                         'odeH2S'    'horizontal flux 2nd order ode S'           'psu.s^-1'
       X2S                         'odeX2S'    'advection flux 2nd order ode S'            'psu.s^-1'       
       rk3_S                       'ode3S'     '3rd order ode for S'                       'psu.s^-1'
       V3S                         'odeV3S'    'vertical flux 3rd order ode S'             'psu.s^-1'
       H3S                         'odeH3S'    'horizontal flux 3rd order ode S'           'psu.s^-1'
       X3S                         'odeX3S'    'advection flux 3rd order ode S'            'psu.s^-1'            
       rk4_S                       'ode4S'     '4th order ode for S'                       'psu.s^-1'
       V4S                         'odeV4S'    'vertical flux 4th order ode S'             'psu.s^-1'
       H4S                         'odeH4S'    'horizontal flux 4th order ode S'           'psu.s^-1'
       X4S                         'odeX4S'    'advection flux 4th order ode S'            'psu.s^-1'             
       };
   
% Write data to file
    
[Sdiag.file, Sdiag.ncid, Sdiag.vid] = diagoutput(Grd, Sdiag, it, datatodiag, tidx(1));

if it == tidx(1)
   cu = onCleanup(@() closefiles([Sdiag.ncid]));
end
   
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


