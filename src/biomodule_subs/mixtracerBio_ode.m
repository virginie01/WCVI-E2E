function [db, Flx] = mixtracerBio_ode(time, bio, In, Grd, ismixed, mld, ent, lb, Sig, tauy, xfil, buoy, names)
% MIXTRACERBIO_ODE Computes mixing and advection for biological tracers in WCVIE2E
%
% This function calculates vertical mixing, horizontal mixing, and lateral
% advection of biological tracers at a given simulation time step.
%
% INPUTS:
%
% time     - Scalar. Time elapsed since simulation start (in seconds)
% bio      - [nz x nx x nbsv] Biological tracer concentrations (mol N m⁻³)
% In      :  structure holding user-supplied input variables
% Grd      - Struct. Grid and temporal data
% ismixed  - [nbsv x 1] Logical vector, true for planktonic (mixed) groups
% mld      - Struct. Mixed layer depth (m) [time x nx]
% ent      - Struct. Entrainment rate (s⁻¹) [time x nx]
% lb       - Nested Struct. Lateral boundary concentrations for each group. 
%            First level = "zooplanktonic" group (i.e. mixed). 
%            Second level:
%                 t:    Grd.time(it):In.datadt:Grd.time(it+1)
%                 o:    1 x 14 array specifying the boundary
%                 data: length(t) x length(o) array holding the data for the current time steps
%                 Data: Column 1 = open ocean upper layer; Column 2 = open ocean lower layer; 
%                 Column 3 = rain shelf; Column 4: rain slope; Column 5 = freshwater from run-offs; 
%                 Column 6 = VICC shelf ul; Column 7 = VICC shelf ll; Column 8 = sbc shelf ul; 
%                 Column 9 = sbc shelf ll'; Column 10 = sbc slope ul; Column 11 = dc shelf ul; 
%                 Column 12 = dc shelf ll; Column 13 = dc slope ul; Column 14 = cu slope ll
% Sig      - [nz x nx] Density field (kg m⁻³)
% tauy     - Struct. North-south wind stress (N m⁻²) [time x nx]
% xfil     - Struct. Cross-shelf transport filter (see advection.m)
% buoy     - Struct. Buoyancy and alongshore currents (see advection.m)
% names    - Cell array [nbsv x 3]. Column 1: short name, 2: long name, 3: units
%
% OUTPUTS:
% db  - [nz x nx x nbsv] Rate of change of tracer concentrations (mol N m⁻³ s⁻¹)
% Flx - Struct with tracer flux fields:
%          V, H, X, CS, R, P, VICC, DC, SBC, CU – all [nbsv+2 x nbsv+2 x nz x nx]
%          and mol N m⁻³ s⁻¹. Here sources and sinks are the same variables
%
% Copyright (c) 2026 Virginie Bornarel
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.


%% Initialize
[nz, nx, nbsv] = size(bio);

% Preallocate flux matrices
[V, H, X, CS, R, P, VICC, DC, SBC, CU] = deal(zeros(nbsv+2, nbsv+2, nz, nx));

% Pad ismixed to include 2 additional diagnostic groups (e.g., non-biological)
ismixed = [ismixed; 0; 0];
idx = find(ismixed);
n = length(idx);

%% Compute Fluxes
for i=1:n

j = idx(i);
name = names{j,2};

v = verticalmixing(bio(:,:,j), Grd, mld, ent, In, time);
V(j,j,:,:) = reshape(v, [1 1 Grd.nz Grd.nx]);

h = horizontalmixing(bio(:,:,j), Grd, mld, lb.(name), In, time);
H(j,j,:,:) = reshape(h, [1 1 Grd.nz Grd.nx]);

[x, cs, r, p, vicc, dc, sbc, cu] = advection(bio(:,:,j), lb.(name), In, Grd, mld, tauy, xfil, Sig, buoy, time);
X(j,j,:,:) = reshape(x, [1 1 Grd.nz Grd.nx]);
CS(j,j,:,:) = reshape(cs, [1 1 Grd.nz Grd.nx]);
R(j,j,:,:) = reshape(r, [1 1 Grd.nz Grd.nx]);
P(j,j,:,:) = reshape(p, [1 1 Grd.nz Grd.nx]);
VICC(j,j,:,:) = reshape(vicc, [1 1 Grd.nz Grd.nx]);
DC(j,j,:,:) = reshape(dc, [1 1 Grd.nz Grd.nx]);
SBC(j,j,:,:) = reshape(sbc, [1 1 Grd.nz Grd.nx]);
CU(j,j,:,:) = reshape(cu, [1 1 Grd.nz Grd.nx]);

end


%% Creating structure Flx

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

%% Total flux (vertical + horizontal + advection)

fluxtot = Flx.V + Flx.H + Flx.X; % nbsv+2 x nbsv+2 x nz x nx molN.m-3.s-1

%% Derivative of tracer concentration (dC/dt)
db = zeros(nz, nx, nbsv);

for i = 1:n
    
    j = idx(i);    
    db(:,:,j) = reshape(fluxtot(j,j,:,:), [Grd.nz Grd.nx 1 1]);
    
end

%% NaN check
if any(isnan(db(:)))
    warning('mixtracerBio_ode:nanInDbdt', 'NaN in dB/dt');
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

idx=mld.t==time & ent.t == time;


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
%       time  : Scalar. Time elapsed since start of the simulation (s).
%
%   Output:
%
%       H   :  [nz x nx]. Horizontal tracer fluxes (tracer/s)

% -----Physical constant --------

Mh=20./86400;% horizontal mixing coefficient (m/s)

% ----- Initialize --------

H=zeros(Grd.nz, Grd.nx);

% ----- Get time index --------

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
function [X, CS, R, P, VICC, DC, SBC, CU] = advection(tracer, lb, In, Grd, mld, tauy, xfil, ~, buoy, time)
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

% Extract values
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
SBC(1,1)=sbc.data(idx).*(lb.data(idx,7)-tracer(1,1));
DC(1,1)=dc.data(idx).*(lb.data(idx,10)-tracer(1,1));

%shelf LL
CS(2,1)=(((afl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(2,2)-tracer(2,1)))...
    +(((dfl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(1,1)-tracer(2,1)));
VICC(2,1)=(VICCshll./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(lb.data(idx,6)-tracer(2,1));
SBC(2,1)=sbc.data(idx).*(lb.data(idx,8)-tracer(2,1));
DC(2,1)=dc.data(idx).*(lb.data(idx,11)-tracer(2,1));

%slope UL
CS(1,2)=(((afl.*(1-0.514932659196421)/In.dx(1,2))/mld.data(idx,2))*(tracer(1,1)-tracer(1,2)))...
    +(((dfl/In.dx(1,2))/mld.data(idx,2))*(lb.data(idx,1)-tracer(1,2)));
P(1,2)=(prate.data(idx,2)./mld.data(idx,2)).*(lb.data(idx,4)-tracer(1,2));
SBC(1,2)=sbc.data(idx).*(lb.data(idx,9)-tracer(1,2));
DC(1,2)=dc.data(idx).*(lb.data(idx,12)-tracer(1,2));

%slope LL
CS(2,2)=(((afl/In.dx(2,2))/100)*(lb.data(idx,2)-tracer(2,2)));%...
   % +(((dfl/In.dx(2,2))/100)*(tracer(2,1)-tracer(2,2)));
CU(2,2)=cu.data(idx).*(lb.data(idx,13)-tracer(2,2));


% X
X= CS + P + R + VICC + DC + SBC + CU;