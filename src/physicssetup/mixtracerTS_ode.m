function [db, Flx] = mixtracerTS_ode(time, tracer, In, Grd, mld, ...
    ent, lb, Sig, tauy, xfil, buoy, varargin)
%   MIXTRACERTS_ODE - Calculates mixing and advection of a tracer in WCVIE2E_physicalmodel
%   
%   Syntax:
%   [db, Flx] = mixtracerTS_ode(time, tracer, In, Grd, mld, ent, lb, Sig, tauy, xfil, buoy, varargin)
%
%   Inputs:
%
%      time    :  Scalar. elapsed simulation time (s).
%
%      tracer  :  [nz x nx] array. Tracer concentrations at the
%                 current time step (units vary based on tracer).
%
%      In      :  struct. User-supplied input variables
%
%      Grd     :  Struct. Grid and temporal metadata
%
%      mld     :  Struct. Contains mixed layer depth data (m) [ time x nx]
%                 at current time step
%
%      ent     :  struct. Contains entrainment rate data (s^-1) [time x nx]
%                 at current time step
%
%      lb      :  Struct. Lateral boundary concentrations of tracer at 
%                 current time step. Contains data [time x 14]
%                 Data: Column 1 = open ocean upper layer; 
%                 Column 2 = open ocean lower layer; Column 3 = rain shelf; 
%                 Column 4: rain slope; Column 5 = freshwater from run-offs; 
%                 Column 6 = VICC ul sh, Column 7 = VICC ll sh; Column 8 =
%                 SBC shelf ul; Column 9 = SBC shelf ll; Column 10 = SBC slope ul; 
%                 Column 11 = DC shelf ul; Column 12 = DC shelf ll; 
%                 Column 13 = DC slope ul; Column 14 = CU slope ll
%
%      Sig     :  [nz x nx] array. Density (kg.m^-3) at current time step
%
%      tauy    :  Struct. Contains S-N wind stress data (N.m^-2) at current 
%                 time steps [time x nx]
%
%      xfil    : Struct. Upwelling index / total forcing containing data 
%                [time x 1]
%
%      buoy    : Struct. Buoyancy-related fluxes or alongshore currents.
%                Data are [time x nx] or [time x 1)
%
%      varargin: Optional name-value pairs for 'sflux' and 'source'.
%                'sflux' refers to heat fluxes across the surface (W.m-2;
%                [1 x nx]). 'source' refers to heat fluxes throughout the
%                water column (W.m-2; [nz x nx])
%
%    Outputs:
%
%      db    :   [nz x nx] array. Total tracer rate of change (tracer
%                units/s)
%      Flx   :   Struct. All individual mixing/advection contributions

%% Core physical process computation

V = verticalmixing(tracer, Grd, mld, ent, In, time);

H = horizontalmixing(tracer, Grd, mld, lb, In, time);

[X, CS, R, P, VICC, DC, SBC, CU] = advection(tracer, lb, In, Grd, mld, tauy, xfil, buoy, time);


%% Check for consistency in dimensions

if ~isequal(size(tracer), size(V), size(H), size (X), size (CS), size (R),...
         size (P), size (VICC), size (DC), size (SBC), size (CU))
    error('flux coefficient array must be same size as tracer array');
end

%% Structure Flx

Flx.V = V;
Flx.H = H;
Flx.X = X;
Flx.CS = CS;
Flx.R = R;
Flx.P = P;
Flx.VICC = VICC;
Flx.DC = DC;
Flx.SBC = SBC;
Flx.CU = CU;

%% Combine fluxes
db = V+H+X;

%% Optional surface and solar heat fluxes 

nArgs=length(varargin);
if nArgs >0
    %surface heat fluxes
    opt = varargin(cellfun(@ischar,varargin));
    if any(ismember(opt, {'sflux'}))

        %recall structures        
        Ht.Qi=varargin{2};
        Ht.airtmp=varargin{3};
        Ht.dewptT=varargin{4};
        Wnd.Wspd10=varargin{5};
        Ht.Qo=varargin{6};
        
        %call the right time index
        idx = Ht.Qi.t==time & Ht.airtmp.t==time & Ht.dewptT.t==time & ...
            Wnd.Wspd10.t==time & Ht.Qo.t==time & mld.t==time;
        %data that correspond to the current time index
        qi = Ht.Qi.data(idx,:);
        airtmp = Ht.airtmp.data(idx,:);
        dewptT= Ht.dewptT.data(idx,:);
        Wspd10= Wnd.Wspd10.data(idx,:);
        Qo=Ht.Qo.data(idx);
        mlayer=mld.data(idx,:);
        
    % Calculate incident, sensible, latent, and longwave heat fluxes, based
    % on heat forcing/stratification at the start of each time step
        [~, qs, ql, qlw] = calcheat(qi, airtmp, dewptT, tracer(1,:),...
            Wspd10, Qo, Grd.nx);
    % Sensible, latent, and longwave fluxes are applied at the surface
    
    srfhflx = qs + ql + qlw;
    
    % Heat fluxes (W/m^2) are converted to temperature fluxes (deg C/s)
    
    Cp = 4125; % Typical specific heat for water (J kg^-1 degC^-1) @10 deg., 32 ppt)
    srf_Tflx = srfhflx./(Sig(1,:).*Cp.* mlayer);
    
    % first optional output 
    
    Flx.srf_Tflx = [srf_Tflx;0 0;0 0];
    
    % calculate new ode including surface heat fluxes
    db= db+ Flx.srf_Tflx;
    
    end
    
    % solar heat fluxes
    
    opt = varargin(cellfun(@ischar,varargin));
    if any(ismember(opt, {'source'}))
        
        %recall structures        
        Ht.Qi=varargin{8};
        
        %call the right time index
        idx = Ht.Qi.t==time & mld.t==time;
        %data that correspond to the current time index
        qi = Ht.Qi.data(idx,:);
        mlayer=mld.data(idx,:);
        
        %update the layer depths for the current time step
       
        % zp values have to be negative and in decreasing order. Add sign
        % minus
        zp = -[Grd.zp(1,1), Grd.zp(1,2); ...
              mlayer(1), mlayer(2);...
              Grd.zp(3,1), Grd.zp(3,2);...
              Grd.zp(4,1), Grd.zp(4,2)];
          
        dz = [mlayer(1), mlayer(2); ...
              Grd.zp(4,1)-In.dz(3,1)-mlayer(1),Grd.zp(4,2)-In.dz(3,2)-mlayer(2);...
              In.dz(3,1), In.dz(3,2)];
         
        
        % Incident heat flux is distributed throughout the water column. prad
        % is the fraction that is photosynthetically-active, attenuated
        % according to krad1.  The remaining fraction is attenuated according
        % to krad2.  Qi is also adjusted downward based on albedo.
    
        solhflx = zeros(Grd.nz,Grd.nx);
        for i = 1:Grd.nx
        solhflx(:,i) = -(1-In.alb) .* qi(i) * diff(...
        In.prad1*exp(In.krad1*zp(:,i)) + (1-In.prad1)*exp(In.krad2*zp(:,i)));        
        end
        
    % Heat fluxes (W/m^2) are converted to temperature fluxes (deg C/s)
    
    Cp = 4125; % Typical specific heat for water (J kg^-1 degC^-1) @10 deg., 32 ppt)
    
    sol_Tflx = solhflx./(Sig.* Cp.* dz);
    
    % second optional output
    
    Flx.sol_Tflx = sol_Tflx;
    
    % new tracer concentration
    db = db + Flx.sol_Tflx;

    end   
        
end

%% Function VERTICALMIXING
function V = verticalmixing(tracer, Grd, mld, ent, In, time)
% VERTICALMIXING Calculates vertical mixing of a tracer
%
% V = verticalmixing(tracer, Grd, mld, ent, In, time)
%
% This function is derived from Ianson and Allen. 2002. A two-dimensional
% nitrogen and carbon flux model in a coastal upwelling region. 
% Vertical mixing (tracer.s^-1) is estimated from the mixed layer depth, 
% entrainment rates and other physical parameters.
%
% Inputs:
%
%       tracer:  [nz x nx]. Tracer concentrations at the current time step
%
%       Grd   :  Grid structure
%
%       mld   :  Structure with .t and .data (mixed layer depths, m)
%
%       ent   :  Structure with .t and .data (entrainment rates; s-1)
%
%       In    : Structure holding user-supplied input variables
%
%       time   :  Scalar. Time elapsed since start of the simulation (s).
%
%   Output:
%
%       V   :  [nz x nx]. Vertical flux of tracer (tracer/s) at current
%              time step

% ----Physical constants ---------

Mv=0.1./86400;   % vertical mixing velocity [m/s]
dm=4;            % mixing depth below MLD [m]
Hshpp=43;        % pycnocline shelf [m]
Hslpp=73;        % pycnocline slope [m]

% ---- Initialize ---------------

V=zeros(Grd.nz, Grd.nx);

% finding time index corrsponding to current time

idx=mld.t==time & ent.t == time;

% dzi and dzo

dzi = Hshpp - mld.data(idx,1);
if dzi <= dm
    dzi = dm;
end

dzo = Hslpp - mld.data(idx,2);
if dzo <= dm
    dzo = dm;
end


%shelf & upper layer (i.e. mixed layer)
V(1,1)=(((Mv.*dm)./(mld.data(idx,1).*dzi))...
    +max(ent.data(idx,1)./(Hshpp+mld.data(idx,1)),0)).*(tracer(2,1)-tracer(1,1));

%shelf & lower layer
V(2,1)=(((Mv.*dm)./((Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)).*dzi))...
    - min(ent.data(idx,1)./((2.*Grd.zp(4,1))-Hshpp+mld.data(idx,1)),0)).*(tracer(1,1)-tracer(2,1))...
    +(Mv./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(tracer(3,1)-tracer(2,1));

%shelf & demersal layer
V(3,1)=(Mv./In.dz(3,1)).*(tracer(2,1)-tracer(3,1));

%slope & upper layer (i.e. mixed layer)
V(1,2)=(((Mv.*dm)./(mld.data(idx,2).*dzo))...
    +max(ent.data(idx,2)./(Hslpp+mld.data(idx,2)),0)).*(tracer(2,2)-tracer(1,2));

%slope & lower layer
V(2,2)=(((Mv.*dm)./((Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2)).*dzo))...
    - min(ent.data(idx,2)./((2.*Grd.zp(4,2))-Hslpp+mld.data(idx,2)),0)).*(tracer(1,2)-tracer(2,2))...
    +(Mv./(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2))).*(tracer(3,2)-tracer(2,2));

%slope & demersal layer
V(3,2)=(Mv./In.dz(3,2)).*(tracer(2,2)-tracer(3,2));

%% Function HORIZONTALMIXING
function H = horizontalmixing(tracer, Grd, mld, lb, In, time)
% HORIZONTALMIXING Calculates horizontal mixing of a tracer
%
% H = horizontalmixing(tracer, Grd, mld, lb, In, time)
%
% This function is derived from Ianson and Allen (2002). A two-dimensional
% nitrogen and carbon flux model in a coastal upwelling region. 
%
%
% Inputs:
%
%       tracer:  [nz x nx]. Tracer concentration matrix
%
%       Grd   : Structure with grid data
%
%       mld   : Structure with mixed layer depth (fields: .t, .data)
%
%       lb    : Nested structure holding lateral boundary concentrations of tracer 
%              
%               t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%               o:    1 x 14 array specifying the boundary
%               data: length(t) x length(o) array holding the data for the current time steps
%               Data: Column 1 = open ocean upper layer; Column 2 = open ocean lower layer; 
%               Column 3 = rain shelf; Column 4: rain slope; Column 5 = freshwater from run-offs; 
%               Column 6 = VICC ul sh; Column 7 = VICC ll sh; Column 8 = sbc shelf ul; 
%               Column 9 = sbc shelf ll'; Column 10 = sbc slope ul; Column 11 = dc shelf ul; 
%               Column 12 = dc shelf ll; Column 13 = dc slope ul; Column 14 = cu slope ll
%
%       In    : Structure holding user-supplied input variables
%
%
%       time  : Scalar. Time elapsed since start of the simulation (s).
%
%   Output:
%
%       H   :  [nz x nx]. Horizontal tracer fluxes (tracer/s)


% -----Physical constant --------

Mh=20./86400; % horizontal mixing coefficient (m/s)

% ----- Initialize --------

H=zeros(Grd.nz, Grd.nx);

% ----- Get time index --------

idx = mld.t==time & lb.t==time;

% ----- Shelf calculations (x = 1) --------

H(1,1)=(Mh./In.dx(1,1)).*(tracer(1,2)-tracer(1,1));

%shelf & lower layer
H(2,1)=(Mh./In.dx(2,1)).*(tracer(2,2)-tracer(2,1));

%shelf & demersal layer
H(3,1)=(Mh./In.dx(3,1)).*(tracer(2,2)-tracer(3,1));

% ----- Slope calculations (x = 2) --------

%slope & upper layer (i.e. mixed layer)
H(1,2)=((Mh./In.dx(1,2)).*(mld.data(idx,1)./mld.data(idx,2)).*(tracer(1,1)-tracer(1,2)))...
    +((Mh./In.dx(1,2)).*(lb.data(idx,1)-tracer(1,2)));

%slope & lower layer
H(2,2)=((Mh./In.dx(2,2)).*((Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))./(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2))).*(tracer(2,1)-tracer(2,2)))...
    +((Mh./In.dx(2,2)).*(In.dz(3,1)./(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2))).*(tracer(3,1)-tracer(2,2)))...
    +((Mh./In.dx(2,2)).*(lb.data(idx,2)-tracer(2,2)));

%slope & demersal layer
H(3,2)=(Mh./In.dx(3,2)).*(lb.data(idx,2)-tracer(3,2));

%% Function ADVECTION
function [X, CS, R, P, VICC, DC, SBC, CU] = advection(tracer, lb, In, Grd, mld, tauy, xfil, buoy, time)
% ADVECTION Calculates advection of a tracer
%
% X = advection(tracer, lb, In, Grd, mld, tauy, xfil, ~, buoy, time)
%
% Based on: Ianson and Allen. 2002. A two-dimensional nitrogen and carbon 
% flux model in a coastal upwelling region. 
%
% Inputs:
%
%      tracer:  [nz x nx]. Tracer concentrations at current time step
%
%      lb    :  Nested structure holding lateral boundary concentrations of tracer 
%               t:    Grd.time(it):In.datadt:Grd.time(it+1), time
%               o:    1 x 14 array specifying the boundary
%               data: length(t) x length(o) array holding the data for the current time steps
%               Data: Column 1 = open ocean upper layer; Column 2 = open ocean lower layer; 
%               Column 3 = rain shelf; Column 4: rain slope; Column 5 = freshwater from run-offs; 
%               Column 6 = VICC ul sh; Column 7 = VICC ll sh; Column 8 = sbc shelf ul; 
%               Column 9 = sbc shelf ll'; Column 10 = sbc slope ul; Column 11 = dc shelf ul; 
%               Column 12 = dc shelf ll; Column 13 = dc slope ul; Column 14 = cu slope ll
%
%      In    :  Structure holding user-supplied input variables
%
%      Grd   :  Structure holding spatial and temporal grid data
%
%      mld   :  Structure holding mixed layer depth data (m). Data:
%               datant+1 x nx array.
%
%      tauy  :  Structure holding S-N wind stress data (N.m^-2).
%               Data: datant+1 x nx array.
%               datant+1 x 1 array if from tauy2
%
%      xfil  :  Structure holding the total upwelling forcing based on 
%               currents from Debby's file
%               datant+1 x 1 array
%
%      buoy :  Nested structure holding buyoancy fluxes and alongshore
%              currents. The sub-structures are PRATE/SBC/DC/UC. PRATE and 
%              CU expressed in m/s. SBC and DC expressed in s-1
%              t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%              data: [length(t) x 1 or 2] array holding the data for the current time steps 
%
%      time  :  Scalar. Time elapsed since start of the simulation (s).
%
%   Output variable:
%
%       X     :  nz x nx array. Flux(es) of tracer due to advection
%                that feed each spatial box (tracer.s^-1) at the current
%                time step
%

prate = buoy.PRATE;
sbc = buoy.SBC;
dc = buoy.DC;
cu = buoy.CU;

%% Indexing

idx = lb.t==time & mld.t==time & tauy.t==time & prate.t==time & sbc.t==time &...
dc.t ==time & cu.t ==time & xfil.t==time;


%% Extract values

force = xfil.data(idx);

if force >0
    a = force;
    d=0;
else
    a=0;
    d=-force;
end

afl=a.*In.dx(1,1);
dfl=d.*In.dx(1,1);

%% Initialize flux arrays
CS=zeros(Grd.nz, Grd.nx);R=zeros(Grd.nz, Grd.nx);
P=zeros(Grd.nz, Grd.nx);VICC=zeros(Grd.nz, Grd.nx);DC=zeros(Grd.nz, Grd.nx);
SBC=zeros(Grd.nz, Grd.nx);CU=zeros(Grd.nz, Grd.nx);

%% Runoff, VICC (converted from m/s)
Rshul=0.7.*(0.11.*exp(-0.3.*(5-(2.*cos(2.*pi.*((time./86400)+20)./365))))-...
    0.06.*exp(-0.6.*(10-(6.*cos(2.*pi./365.*((time./86400)-150))))));%model modified
Rshul=Rshul./86400;

VICCshul=0.5.*exp(-0.6.*(5-(2.*cos(2.*pi./365.*((time./86400)-150)))));%model
VICCshul=VICCshul./86400;

VICCshll = VICCshul.*(50-mld.data(idx,1))./mld.data(idx,1);


%% Shelf Upper Layer
CS(1,1)=(((afl/In.dx(1,1))/mld.data(idx,1))*(tracer(2,1)-tracer(1,1)))...
    +(((dfl/In.dx(1,1))/mld.data(idx,1))*(tracer(1,2)-tracer(1,1)));
P(1,1)=(prate.data(idx,1)./mld.data(idx,1)).*(lb.data(idx,3)-tracer(1,1));
R(1,1)=(Rshul./mld.data(idx,1)).*(lb.data(idx,5)-tracer(1,1));
VICC(1,1)=(VICCshul./mld.data(idx,1)).*(lb.data(idx,6)-tracer(1,1));
SBC(1,1)=sbc.data(idx).*(lb.data(idx,8)-tracer(1,1));
DC(1,1)=dc.data(idx).*(lb.data(idx,11)-tracer(1,1));

%% Shelf Lower Layer
CS(2,1)=(((afl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(2,2)-tracer(2,1)))...
    +(((dfl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(1,1)-tracer(2,1)));
VICC(2,1)=(VICCshll./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(lb.data(idx,7)-tracer(2,1));
SBC(2,1)=sbc.data(idx).*(lb.data(idx,9)-tracer(2,1));
DC(2,1)=dc.data(idx).*(lb.data(idx,12)-tracer(2,1));

%% Slope Upper Layer
CS(1,2)=(((afl.*(1-0.514932659196421)/In.dx(1,2))/mld.data(idx,2))*(tracer(1,1)-tracer(1,2)))...
    +(((dfl/In.dx(1,2))/mld.data(idx,2))*(lb.data(idx,1)-tracer(1,2)));
P(1,2)=(prate.data(idx,2)./mld.data(idx,2)).*(lb.data(idx,4)-tracer(1,2));
SBC(1,2)=sbc.data(idx).*(lb.data(idx,10)-tracer(1,2));
DC(1,2)=dc.data(idx).*(lb.data(idx,13)-tracer(1,2));

%% Slope Lower Layer
CS(2,2)=(((afl/In.dx(2,2))/100)*(lb.data(idx,2)-tracer(2,2)));%...
   % +(((dfl/In.dx(2,2))/100)*(tracer(2,1)-tracer(2,2)));
CU(2,2)=cu.data(idx).*(lb.data(idx,14)-tracer(2,2));


%% Total Advection Flux
X= CS + P + R + VICC + DC + SBC + CU;