function varargout = physicalmodel(varargin)
%PHYSICALMODEL Run a 2-D physical–biological box model for a coastal upwelling system
%
%   physicalmodel(outputfolder)
%   physicalmodel(outputfolder, param1, val1, ...)
%
% DESCRIPTION
%   This function is the main driver for a 2-D box model representative of a
%   coastal upwelling system, based on Ianson & Allen (2002). It integrates:
%
%     - Physical tracers (temperature, salinity)
%     - Optional biological / E2E modules
%     - Vertical mixing and advection
%     - Archiving of model outputs
%
%   The model can be run in:
%     (1) physics-only mode
%     (2) coupled physical–biological mode
%
% INPUTS
%
%   REQUIRED:
%   ---------
%
%   outputfile:   string, base name of output folder for simulation(s).  
%                 The output folder will contain a table with all model
%                 outputs, including all state variables (see
%                 WCVI-E2E_run.m and associated documentation)
%
%   MODEL GRID:
%   -----------
%
%   dz:         A 2D nz*nx array, the thickness of each box in the vertical
%               dimension
%   
%   dx:         a 2D nz*nx array, the width of each box in the horizontal
%               dimension        
% 
%   dt:         the model time step (seconds) [10800 ~ 3 hrs]
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
%               force calculations) [49.5]
%
%   whgt:       elevation above sea level where wind forcing data was
%               measured (m). [10]  
%
%   INITIAL CONDITIONS:
%   -------------------
%        
%   t_input:        Initial temperature profile for simulation.
%                   Data is an nz x nx matrix with initial temperatures.
%              
%   s_input:        Initial salinity profile for simulation.
%                   Data is an nz x nx matrix with initial salinities.
%
%   EXTERNAL FORCING:
%   -----------------
%   See documentation in forcingdata.m for a description of all external
%   forcing datasets.
%   External forcing data are loaded and imported at each simulation time
%   step through the forcingdata.m to avoid memory issues.
%   
%   ARCHIVING:
%   ----------
%
%   tarch:          the archiving interval. Data is averaged over the
%                   specified tarch interval. Data is averaged over tarch
%                   seconds. [3600 x 24 x 30]
%
%   beginarchive:   datenumber indicating what model date to begin recording
%                   to an output file, or NaN to indicate that archiving begins
%                   immediately. [NaN]
%
%   endarchive:     datenumber indicating what model date to stop recording
%                   output to an output file, or NaN to indicate archiving until
%                   the simulation ends. [NaN]
%
%   stopafterinit: logical scalar.  If true, simulation is terminated after
%                  the initialization process, and no forward integration is
%                  performed.  This is for debugging purposes.  Using this
%                  option also leads to different outputs if physicalmodel is
%                  called with an output variable (see below) [false]
%   BIOLOGY:
%   --------
%
%   biofun:         Function handle to biological/e2e module.  If empty, 
%                   the model will run without any biology. [empty] 
%                   If the biological module 'biomodel' is called, see
%                   parsebioinput.m for required biological datasets, IC,
%                   and parameters
%
%   OTHER
%   -----
%
%   verbose:        logical scalar.  If true, progress statements will be
%                   printed to the screen.  If false, nothing will be printed.
%                   [true]
%
% OUTPUTS
%   If an output argument is requested:
%
%     physicalmodel(...) -> returns Arch.file
%
%   If 'stopafterinit' is true:
%
%     physicalmodel(...) -> returns a struct containing initialized state
%
%

% ------------------------------------------------------------------------
%  Setup and input parsing
%  ------------------------------------------------------------------------
mltimer = tic;

In = parseinput(varargin{:});

% Handle early-exit (initialization-only) mode
if isfield(In, 'stopafterinit') % flag to stop after initialization
    stopafterinit = true;
    In = rmfield(In, 'stopafterinit');
else
    stopafterinit = false;
end

if In.verbose
    filename_no_ext = erase(In.outputfile, '.nc');
    fprintf('\n--------------------\n2D Physical Model\n--------------------\n\n');
    fprintf('%s\n\n', filename_no_ext);  
end

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------

if In.verbose
    fprintf('Initializing model...\n');
end

[Grd, TS, Arch] = initialize(In);

% Initialize biology or physics-only diagnostics
if In.hasbio
    [Bio.bioin, Bio.bio, Bio.ismixed, Bio.Params, Bio.names, ...
        Bio.diagnames] = feval(In.biofun, 'init', In, Grd);
    Bio.ndiag = size(Bio.diagnames,1);
    if Bio.ndiag > 0
        Bio.diag = zeros(Grd.nz,Grd.nx,Bio.ndiag); % placeholder KAK note: Charlie's new version runs biofun here to get ititial conditions
    end  
else
    [Phy.diagnames, Phy.flux] = phydiagnostic();
    Phy.ndiag = size(Phy.diagnames,1);
    Phy.diag = zeros(Grd.nz,Grd.nx,Phy.ndiag);
end

% Optional early stop for debugging
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

% Set up progress printing strings
simdatestr = datestr(datenum(Grd.start_date) + Grd.time./86400, '  mm/dd/yyyy\n');
dummystr = [repmat(' ', 1, size(simdatestr,2))-1 '\n'];
erasestr = repmat('\b', 1, size(simdatestr,2));

tidx = 1:Grd.nt;


% -------------------------------------------------------------------------
%  Time stepping loop
%  ------------------------------------------------------------------------
for it=tidx
    
%% Load forcing variables for the current time step
if In.hasbio
[Ht, Wnd, VMxng, LB, Buoy, Bio] = forcingdata(Grd, it, In, Bio);
else
[Ht, Wnd, VMxng, LB, Buoy] = forcingdata(Grd, it, In);
end    

%% set parameter structures for physics and grid for the current time step
PhysParams = TS;
PhysParams.par24=Ht.QI.data(Ht.QI.t == Grd.time(it),:).*In.prad1;%1 x nx array
PhysParams.kpar = In.krad1;

mld = VMxng.MLD.data(VMxng.MLD.t == Grd.time(it),:);
dz =  [mld(1),mld(2); ...
      Grd.zp(4,1)-In.dz(3,1)-mld(1),Grd.zp(4,2)-In.dz(3,2)-mld(2);...
      In.dz(3,1), In.dz(3,2)];
zp = [0 0; cumsum(dz,1)];
z = (zp(1:end-1,:) + zp(2:end,:))./2;
   
if In.hasbio
GrdParams = struct('z', z, 'dz', dz, 'x', Grd.x, 'dx', In.dx, 'nz', Grd.nz, 'nx', Grd.nx, 'area',...
    Bio.bioin.area, 't', Grd.time(it), 'dt', In.dt, 'sdate', Grd.start_date);
else
 GrdParams = struct('z', z, 'dz', dz, 'x', Grd.x, 'dx', In.dx, 'nz', Grd.nz, 'nx', Grd.nx,...
     't', Grd.time(it), 'dt', In.dt, 'sdate', Grd.start_date);
end
   

%% Archiving process for the current time step

% Base physical variables to archive
datatoarchive = {...
 TS.T            'temp'      'temperature'       'deg C'
 TS.S            'sal'       'salinity'          'psu'
 TS.Sig          'sig'       'density'           'kg m^-3'
 };    

% addtional variables to archive    
if In.hasbio
      biodata = num2cell(cat(3,Bio.bio,Bio.diag),[1 2]); %1 x 1 x (nbio + ndiag) cell array. Each cell contains a 3-by-2 array of Bio.bio or Bio.diag
      biodata = reshape(biodata,size(biodata,3),1,1); % (nbio + ndiag) x 1 x 1 cell array. Each cell contains a 3-by-2 array of Bio.bio or Bio.diag
      datatoarchive = [datatoarchive; [biodata [Bio.names; Bio.diagnames]]];
else
    diagdata = num2cell(Phy.diag, [1 2]);
    diagdata = reshape(diagdata,size(diagdata,3),1,1);
    datatoarchive = [datatoarchive; [diagdata Phy.diagnames]];
end
               
% Archiving
for io = 1:length(Arch)   
     [Arch(io).avg, Arch(io).datatoarchive, Arch(io).file] = archivedata(GrdParams, Arch(io), it, datatoarchive, tidx(1)); 
end

%% Print simulation date to screen   
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

%% Advance all biological state variables by one time step

if In.hasbio

% Diapause boundary-condition updates
if Bio.Params.diapause    
    Bio.LB = diapauseLBbio(Bio.LB, Bio.bio, Bio.names, Bio.Params);    
end
                   
% Horizontal advection / mixing of planktonic variables
if Bio.Params.diapause
    Bio.bio(:,:,Bio.Params.idx.zl) = 0;
end

[Bio.bio, dbmix, Flxmix] = rk4_Bio(@mixtracerBio_ode, it, Bio.bio, In, Grd, Bio.ismixed, VMxng.MLD,... 
                                   VMxng.ENT, Bio.LB, TS.Sig, Wnd.TAUY, Wnd.XFIL, Buoy, Bio.names);

if Bio.Params.diapause
    fx = fieldnames(Flxmix);
    for ii = 1:length(fx)
    Flxmix.(fx{ii})(Bio.Params.idx.zl,Bio.Params.idx.zl,:,:) = Flxmix.(fx{ii})(Bio.Params.idx.zl1,Bio.Params.idx.zl1,:,:) + Flxmix.(fx{ii})(Bio.Params.idx.zl2,Bio.Params.idx.zl2,:,:);
    end
end
        
if any(isnan(Bio.bio(:)))
   warning('ML: NaN in biology')
end


% vertical movement (sinking/swimming)
Bio.wsink = feval(In.biofun, 'vertmove', Bio.Params, GrdParams);            

[Bio.bio, dbvflx, Flxvflx] = verticalflux(Bio.bio, Bio.wsink, In.dt, GrdParams);
             
if Bio.Params.diapause
    fx = fieldnames(Flxvflx);
    for ii = 1:length(fx)
    Flxvflx.(fx{ii})(Bio.Params.idx.zl,Bio.Params.idx.zl,:,:) = Flxvflx.(fx{ii})(Bio.Params.idx.zl1,Bio.Params.idx.zl1,:,:) + Flxvflx.(fx{ii})(Bio.Params.idx.zl2,Bio.Params.idx.zl2,:,:);
    end
end

if any(isnan(Bio.bio(:)))
   warning('ML: NaN in biology')
end

% Biological source/sink terms
if Bio.bioin.isnem
  [Bio.bio, dbbio, Flxbio, Diag] = feval(In.biofun, 'sourcesink', Bio.bioin.isnem, Bio.bio,...
      PhysParams, Bio.Params, GrdParams, Bio.O2, Arch, it);                                   
else
  [Bio.bio, dbbio, Flxbio, Diag] = feval(In.biofun, 'sourcesink', Bio.bioin.isnem, Bio.bio,...
      PhysParams, Bio.Params, GrdParams, Bio.O2, Arch, it, Bio.BIOTS);
end  
                                     
if any(isnan(Bio.bio(:)))
   warning('ML: NaN in biology')
end

% Sum tendencies and merge fluxes 
db = dbmix + dbvflx + dbbio;
Flx = mergestruct(Flxbio, Flxmix, Flxvflx);

% Construct Bio.diag (db + selected fluxes + diagnostics)
nd = size(Bio.Params.flux,1);

fluxes = zeros(size(z,1), size(z,2), nd); % nz x nx x nbfluxes
for id = 1:nd   
    fluxes(:,:,id) = Flx.(Bio.Params.flux{id,1})(Bio.Params.flux{id,2},Bio.Params.flux{id,3},:,:); % mol N.m-3.s-1 nz x nx x nbfluxes same order as in listfluxes
end

ndiag = 15;
if isempty(Diag)
    diag = zeros(size(z,1), size(z,2), ndiag);
else
    diag = struct2cell(Diag);
    diag = cat(3, diag{:});%nx x nx x ndiag
end

Bio.diag = cat(3, db, fluxes, diag);

end

%% Advance T, S, Sig state variables by one time step

% Temperature (T)
[TS.T, dbT, FlxT] = rk4_TS(@mixtracerTS_ode, it, ...
TS.T, In, Grd, VMxng.MLD, VMxng.ENT, LB.TMP, TS.Sig, Wnd.TAUY, Wnd.XFILPHY, Buoy, ...
'sflux', Ht.QI, Ht.AIRTMP, Ht.DEWPTT, Wnd.WSPD10, Ht.QO,'source',Ht.QI);

% Salinity (S)
[TS.S, dbS, FlxS] = rk4_TS(@mixtracerTS_ode, it, TS.S, In, Grd, ...
    VMxng.MLD, VMxng.ENT, LB.SAL, TS.Sig, Wnd.TAUY, Wnd.XFILPHY, Buoy);

% Update density
TS.Sig = sw_dens0(TS.S, TS.T);

% Physics-only diagnostics
if ~In.hasbio
    db = cat(3,dbT,dbS);
    fld = fieldnames(FlxT);
    for i = 1:length(fld)
        Flx.(fld{i})(1,1,:,:)=reshape(FlxT.(fld{i}),[1 1 Grd.nz Grd.nx]);
        if ~strcmp(fld{i},'sol_Tflx')&&~strcmp(fld{i},'srf_Tflx')
            Flx.(fld{i})(2,2,:,:)=reshape(FlxS.(fld{i}),[1 1 Grd.nz Grd.nx]);
        end
    end

    nd = size(Phy.flux,1);
    fluxes = zeros(size(z,1),size(z,2),nd);
    for id = 1:nd
        fluxes(:,:,id)=Flx.(Phy.flux{id,1})(Phy.flux{id,2},Phy.flux{id,3},:,:);
    end
Phy.diag = cat(3,db,fluxes);
end

end

%% Finalization
runtime = toc(mltimer);

% Print loop runtime
if In.verbose
	fprintf('\nMixed layer model completed successfully: %f s\n', runtime);
end

% Return archived output if requested
if nargout == 1
    varargout{1} = Arch.file;
end

% Cleanup
if In.hasbio && isfield(Bio.Params, 'cfid')
    fclose(Bio.Params.cfid);
end

function closefiles(id)
for ii = 1:length(id)
    netcdf.close(id(ii));
end


