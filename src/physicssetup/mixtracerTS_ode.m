function [db, Flx] = mixtracerTS_ode(time, tracer, In, Grd, mld, ...
    ent, lb, Sig, tauy, xfil, buoy, varargin)
%   MIXTRACER Calculates mixing of a tracer in WCVIE2E_physicalmodel
%
%   dxdt = mixtracer_ode(tracer, In, Grd, mld, ent, lb, Sig, tauy,prate, it, varargin)
%
%   Input variables:
%
%      time    :  Scalar. Time elapsed since start of the simulation (s).
%
%      tracer  :  nz x nx array. Tracer concentrations in each spatial box at the
%                 current time step (units vary based on tracer).
%
%      In      :  structure holding user-supplied input variables
%
%      Grd     :  Structure holding spatial and temporal grid data for
%                 WCVIE2E_physicalmodel simulations
%
%      mld     :  Structure holding mixed layer depth data (m) at current time steps
%
%      ent     :  structure holding entrainment rate data (s^-1) at current time steps
%
%      lb      :  Structure holding lateral boundary concentrations of 
%                 tracer (units vary based on tracer) at current time steps
%                 Data: Column 1 = open ocean upper layer; Column 2 = open ocean lower layer; 
%                 Column 3 = rain shelf; Column 4: rain slope; Column 5 = freshwater from run-offs; 
%                 Column 6 = VICC
%
%      Sig     :  nz x nx array holding the density (kg.m^-3) values in each
%                 spatial box at the current time step
%
%      tauy    :  Structure holding S-N wind stress data (N.m^-2) at current time steps
%
%      prate   :  Structure holding precipitation rates at current time steps. Expressed 
%                 in m.s^-1.
%
%
%   Optional input variables passed as parameter + structures needed to
%   calculate the parameter (i.e. relevant for temperature):
%
%      sflux   :   heat fluxes across the surface. Will be calculated as an 1 x nx array. 
%                  Expressed in W.m-2                 
%
%      bflux   :   heat fluxes across the seafloor. Will be calculated as an 1 x nx array. 
%                  Expressed in W.m-2
%
%      source  :   heat fluxes throughout the water column in all boxes. Will be calculated as an 
%                  nz x nx array. Expressed in W.m-2
%
%    Output variable:
%
%      dxdt    :   nz x nx array. ODE for tracer concentrations(tracer unit.s^-1).

%% Calculating each component of mixing and advection

V = verticalmixing(tracer, Grd, mld, ent, In, time);

H = horizontalmixing(tracer, Grd, mld, lb, In, time);

[X, CS, R, P, VICC, DC, SBC, CU] = advection(tracer, lb, In, Grd, mld, tauy, xfil, buoy, time);


%% Check input

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

%% Calculating new tracer concentrations
db = V+H+X;

%% Deal with optional arguments: name of parameter to be calculated 
%% followed by structures needed to calculate the parameter.  

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
% VERTICALMIXING calculates vertical mixing of a tracer in the
% WCVIE2E_physicalmodel
%
% V = verticalmixing(tracer, Grd, mld, ent, In, time)
%
% This function is derived from Ianson and Allen. 2002. A two-dimensional
% nitrogen and carbon flux model in a coastal upwelling region. 
% Vertical mixing (tracer.s^-1) is estimated from the mixed layer depth, 
% entrainment rates and other physical parameters.
%
%   Input variables:
%
%       tracer:  nz x nx array. Tracer concentrations in each spatial box 
%                at the current time step (units vary based on tracer)
%
%       Grd   :  Structure holding spatial and temporal grid data for
%                WCVIE2E_physicalmodel simulations
%
%       mld   :  Structure holding mixed layer depth data (m) at current time steps
%
%       ent   :  Structure holding entrainment rate data (s^-1) at current time steps
%
%       In    :  structure holding user-supplied input variables
%
%       time   :  Scalar. Time elapsed since start of the simulation (s).
%
%   Output variable:
%
%       V   :  nz x nx array. Flux(es) of tracer in the vertical dimension
%              that feed each spatial box (tracer.s^-1)at the current time
%              step.

%Physical parameters

Mv=0.1./86400;   % vertical mixing (m.s^-1)- 0.2 in paper but 0.1 in model
dm=4; %depth of mixing below hu (m)- 4 m in model and 2 in paper-
Hshpp=43; %permanent pycnocline shelf (m)
Hslpp=73; %permanent pycnocline slope (m)


% Creation of the matrix V

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
% HORIZONTALMIXING calculates horizontal mixing of a tracer in the
% WCVIE2E_physicalmodel
%
%H = horizontalmixing(tracer, Grd, mld, lateralforcing, In, time)
%
%This function is derived from Ianson and Allen. 2002. A two-dimensional
%nitrogen and carbon flux model in a coastal upwelling region. 
%
%
%   Input variables:
%
%       tracer:  nz x nx array. Tracer concentrations in each spatial box 
%                at the current time step (units vary based on tracer)
%
%       Grd   : Structure holding spatial and temporal grid data for
%               WCVIE2E_physicalmodel simulations
%
%       mld   : Structure holding mixed layer depth data (m) at current time steps
%
%       lb    : Nested structure holding lateral boundary concentrations of tracer 
%               (units vary based on tracer)
%               t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%               o:    1 x 6 array specifying the boundary
%               data: length(t) x length(o) array holding the data for the current time steps
%               Data: Column 1 = open ocean upper layer; Column 2 = open ocean lower layer; 
%               Column 3 = rain shelf; Column 4: rain slope; Column 5 = freshwater from run-offs; 
%               Column 6 = VICC
%
%       In    : Structure holding user-supplied input variables
%
%       time  : Scalar. Time elapsed since start of the simulation (s).
%
%   Output variable:
%
%       H   :  nz x nx array. Flux(es) of tracer in the horizontal dimension
%              that feed each spatial box (tracer.s^-1) at the current time
%              step

%Physical parameters

Mh=20./86400;% horizontal mixing m.s-1

% Creation of the matrix H

H=zeros(Grd.nz, Grd.nx);

% logical vector indicating current time 

idx = mld.t==time & lb.t==time;

%shelf & upper layer (i.e. mixed layer)
H(1,1)=(Mh./In.dx(1,1)).*(tracer(1,2)-tracer(1,1));

%shelf & lower layer
H(2,1)=(Mh./In.dx(2,1)).*(tracer(2,2)-tracer(2,1));

%shelf & demersal layer
H(3,1)=(Mh./In.dx(3,1)).*(tracer(2,2)-tracer(3,1));

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
% ADVECTION calculates advection of a tracer in the
% WCVIE2E_physicalmodel
%
%   X = advection(tracer, lb, In, Grd, mld, tauy, Sig, prate, it, time)
%
%   This function is derived from Ianson and Allen. 2002. A two-dimensional
%   nitrogen and carbon flux model in a coastal upwelling region. 
%
%   Input variables:
%
%      tracer:  nz x nx array. Tracer concentrations in each spatial box at the
%               current time step (units vary based on tracer)
%
%      lb    :  Nested structure holding lateral boundary concentrations of tracer 
%               (units vary based on tracer)
%               t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%               o:    1 x 13 array specifying the boundary
%               data: length(t) x length(o) array holding the data for the current time steps
%               Data: Column 1 = open ocean upper layer; Column 2 = open ocean lower layer; 
%               Column 3 = rain shelf; Column 4: rain slope; Column 5 = freshwater from run-offs; 
%               Column 6 = VICC; Column 7 = sbc shelf ul; Column 8 = sbc shelf ll';
%               Column 9 = sbc slope ul; Column 10 = dc shelf ul; Column 11
%               = dc shelf ll; Column 12 = dc slope ul; Column 13 = cu
%               clope ll
%
%      In    :  Structure holding user-supplied input variables
%
%      Grd   :  Structure holding spatial and temporal grid data for
%               WCVIE2E_physicalmodel simulations
%
%      mld   :  Structure holding mixed layer depth data (m). Data:
%               datant+1 x nx array.
%
%      tauy  :  Structure holding S-N wind stress data (N.m^-2).
%               Data: datant+1 x nx array.
%               datant+1 x 1 array if from tauy2
%
%      xfil  :  Structure holding the remote upwelling index (m.s-1) or the
%               total upwelling forcing based on currents from Debby's file
%               datant+1 x 1 array
%
%      Sig   :  nz x nx array holding the density (kg.m^-3) values in each
%               spatial box at the current time-step
%
%      buoy :  Nested structure holding buyoancy fluxes and alongshore
%              currents. The sub-structures are PRATE/SBC/DC/UC. PRATE and 
%              CU expressed in m/s. SBC and DC expressed in s-1
%              t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%              o:    1 x 13 array specifying the boundary
%              data: length(t) x length(o) array holding the data for the current time steps 
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
% logical vector indicating current time 

 idx = lb.t==time & mld.t==time & tauy.t==time & prate.t==time & sbc.t==time &...
 dc.t ==time & cu.t ==time & xfil.t==time;

 
 %% If upwelling forcing calculated based on winds: Ekman velocity + Remote upwelling index from Ze
%-----------------------------
%Calculating Ekman velocity
%-----------------------------
%Physical parameters

%cf = 2*7.2921*10.^(-5)*sin(In.Lat*(pi/180)); %Coriolis parameter (rad.s^-1)
%Rd = 20000; %internal Rossby radius of deformation (m)

% calculating local a/d

%uek=tauy.data(idx,1)/(Sig(1,1)*cf*Rd); %uek in m.s^-1

%if uek > 0
%    alocal=0;
%    dlocal=uek; % northward wind (+ ve) causes downwelling (-ve)
%else
%    alocal=-uek;
%    dlocal=0; % southward wind (-ve) causes upwelling (+ve)
%end
%------------------------
% calculating remote a/d
%------------------------
%uremote = xfil.data(idx);

%if uremote >0
%    aremote=0;
%    dremote=uremote;
%else
%    aremote=-uremote;
%    dremote=0;
%end

%a = alocal + aremote;
%d = dlocal + dremote;

%% if upwelling forcing calculated based on currents from Debby's files

force = xfil.data(idx);

if force >0
    a = force;
    d=0;
else
    a=0;
    d=-force;
end

afl=a.*In.dx(1,1);%a*shelf width
dfl=d.*In.dx(1,1);%d*shelf width

%creating CS(cross-shore), R(run-offs), P(precipitation), VICC, 
%DC(Davidson current), SBC(shelf break current), CU(California Undercurrent)
%matrices
CS=zeros(Grd.nz, Grd.nx);R=zeros(Grd.nz, Grd.nx);
P=zeros(Grd.nz, Grd.nx);VICC=zeros(Grd.nz, Grd.nx);DC=zeros(Grd.nz, Grd.nx);
SBC=zeros(Grd.nz, Grd.nx);CU=zeros(Grd.nz, Grd.nx);

%%Buoyancy fluxes

%flux per length of coastline for rainfall- Both publications and model have
% similar cycles even though the formula differ- Formula from model here-

%P= 0.03.*exp(-0.3.*(5-(2.*cos(2.*pi.*((time./86400)+20)/365))))...
%   -0.03.*exp(-0.6.*(10-(6.*cos(2.*pi/365.*((time/86400)-150)))));
%P=P./86400;

%flux per length of coastline for terrigenous runoff: 1)Formula from 
%Debby's model-2)first equation modified to replicate figures found in
%literature (i.e. min runoff = 1 x 10^3 m3/s-->6.8 x 10^-8 m/s;
% max runoff = 5 x 10^3 m3/s -->3.4 x 10^-7 m/s)

%R=0.5.*(0.11.*exp(-0.3.*(5-(2.*cos(2.*pi.*((time./86400)+20)./365))))-...
%    0.11.*exp(-0.6.*(10-(6.*cos(2.*pi./365.*((time./86400)-150)))))); %model
%R=R./86400;

Rshul=0.7.*(0.11.*exp(-0.3.*(5-(2.*cos(2.*pi.*((time./86400)+20)./365))))-...
    0.06.*exp(-0.6.*(10-(6.*cos(2.*pi./365.*((time./86400)-150))))));%model modified
Rshul=Rshul./86400;


%flux per length of coastline for the VICC-1)Formula from 
%Debby's model corresponds to literature 
%(i.e. max runoff = 1.4 x 10^5 m3/s or flux with speed 0.15m/s, current width
% 20km, depth=50m--> gives 1 x 10^-6 m/s in summer in UL)

VICCshul=0.5.*exp(-0.6.*(5-(2.*cos(2.*pi./365.*((time./86400)-150)))));%model
VICCshul=VICCshul./86400;

% to convert VICC in UL to VICC in LL, multiply by (50-hsh,ul)/hsh,ul
VICCshll = VICCshul.*(50-mld.data(idx,1))./mld.data(idx,1);


% X = CS + P + R + VICC + DC + SBC + CU
%shelf UL
CS(1,1)=(((afl/In.dx(1,1))/mld.data(idx,1))*(tracer(2,1)-tracer(1,1)))...
    +(((dfl/In.dx(1,1))/mld.data(idx,1))*(tracer(1,2)-tracer(1,1)));
P(1,1)=(prate.data(idx,1)./mld.data(idx,1)).*(lb.data(idx,3)-tracer(1,1));
R(1,1)=(Rshul./mld.data(idx,1)).*(lb.data(idx,5)-tracer(1,1));
VICC(1,1)=(VICCshul./mld.data(idx,1)).*(lb.data(idx,6)-tracer(1,1));
SBC(1,1)=sbc.data(idx).*(lb.data(idx,8)-tracer(1,1));
DC(1,1)=dc.data(idx).*(lb.data(idx,11)-tracer(1,1));

%shelf LL
CS(2,1)=(((afl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(2,2)-tracer(2,1)))...
    +(((dfl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(1,1)-tracer(2,1)));
VICC(2,1)=(VICCshll./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(lb.data(idx,7)-tracer(2,1));
SBC(2,1)=sbc.data(idx).*(lb.data(idx,9)-tracer(2,1));
DC(2,1)=dc.data(idx).*(lb.data(idx,12)-tracer(2,1));

%slope UL
CS(1,2)=(((afl.*(1-0.514932659196421)/In.dx(1,2))/mld.data(idx,2))*(tracer(1,1)-tracer(1,2)))...
    +(((dfl/In.dx(1,2))/mld.data(idx,2))*(lb.data(idx,1)-tracer(1,2)));
P(1,2)=(prate.data(idx,2)./mld.data(idx,2)).*(lb.data(idx,4)-tracer(1,2));
SBC(1,2)=sbc.data(idx).*(lb.data(idx,10)-tracer(1,2));
DC(1,2)=dc.data(idx).*(lb.data(idx,13)-tracer(1,2));

%slope LL
CS(2,2)=(((afl/In.dx(2,2))/100)*(lb.data(idx,2)-tracer(2,2)));%...
   % +(((dfl/In.dx(2,2))/100)*(tracer(2,1)-tracer(2,2)));
CU(2,2)=cu.data(idx).*(lb.data(idx,14)-tracer(2,2));


% X
X= CS + P + R + VICC + DC + SBC + CU;