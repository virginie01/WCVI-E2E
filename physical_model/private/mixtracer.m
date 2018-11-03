
function newtracer = mixtracer(tracer, In, Grd, mld, ent, lb, wndl, wndr, Sig, it, varargin)
%   MIXTRACER Calculates mixing of a tracer in WCVIE2E_physicalmodel
%
%   newtracer = mixtracer(tracer, In.dt, Grd, mld, ent, In, lb, wndl, wndr, Sig, it, varargin)
%
%   Input variables:
%
%      tracer  :  nz x nx array. Tracer concentrations in each spatial box at the
%                 current time step (units vary based on tracer).
%
%      Grd     :  Structure holding spatial and temporal grid data for
%                 WCVIE2E_physicalmodel simulations
%
%      mld     :  1 x nx array. Mixed layer depth at the current time step (m).
%
%      ent     :  1 x nx array. Entrainment rates at the current time step
%                 (day^-1)
%
%      In      :  structure holding user-supplied input variables
%
%      lb      :  1 x 5 array. Lateral boundary concentrations of tracer 
%                 (units vary based on tracer). Column 1 = open ocean upper
%                 layer; Column 2 = open ocean lower layer; Column 3 = rain;
%                 Column 4 = freshwater from run-offs; Column 5 = VICC
%
%      wndl    :  scalar. Local longshore wind stress (N.m^-2) at the
%                 current time step.
%
%      wndr    :  stucture related to remote wind forcing data
%
%      Sig     :  nz x nx array holding the density (kg.m^-3) values in each
%                 spatial box at the current time step
%
%      it      :  current time step
%
%   Optional input variables (passed as parameter/value pairs):
%
%      sflux   :   1 x nx array. Flux of tracer across the surface interface. 
%                  (tracer unit d^-1)
%
%      bflux   :   1 x nx array. Flux of tracer across the bottom interface.
%                  (tracer unit d^-1)
%
%      sval    :   1 x nx array. Tracer value at the surface, used to force
%                  the surface spatial boxes.(Tracer unit).
%
%      bval    :   1 x nx array. Tracer value along the bottom, used to 
%                  force the bottom spatial boxes. (Tracer unit).
%
%      source  :   nz x nx array, source (or sink) flux of tracer in each
%                  spatial unit(tracer unit d^-1).
%
%      dissipate:  dissipation constant (d^-1)   
%
%
%    Output variable:
%
%      newtracer:  nz x nx array. Tracer concentrations at next time step
%                  (tracer unit).

%% Calculating each component of mixing and advection

V = verticalmixing(tracer, Grd, mld, ent, In);
H = horizontalmixing(tracer, Grd, mld, lb, In);
X = advection(tracer, lb, In, Grd, mld, wndl, wndr, Sig, it);

%% Check input

if ~isequal(size(tracer), size(V), size(H), size (X))
    error('flux coefficient array must be same size as tracer array');
end

%% Defaults for optional parameters

A.sflux     = zeros(1, Grd.nx);
A.bflux     = zeros(1, Grd.nx);
A.sval      = zeros(1, Grd.nx);
A.bval      = zeros(1, Grd.nx);
A.source    = zeros(Grd.nz,Grd.nx);
A.dissipate = zeros(Grd.nz,Grd.nx);

%% Calculating new tracer concentrations
newtracer = tracer + (V+H+X).*In.dt;

%% Adjust new tracer concentrations when source (or sink) flux of tracer

if ~isnan(A.source)
    newtracer = tracer + (V+H+X).*In.dt + A.source.*In.dt;
end
    
%% Adjust concentrations for flux boundary conditions

if ~isnan(A.sflux)
    newtracer(1,:)= tracer(1,:) + (V(1,:)+H(1,:)+X(1,:)).*In.dt + A.sflux.*In.dt;
end

if ~isnan(A.bflux)
    newtracer(max(Grd.nz),:)=tracer(max(Grd.nz),:) + (V(max(Grd.nz),:)+H(max(Grd.nz),:)+X(max(Grd.nz),:)).*In.dt - A.bflux.*In.dt;
end

%% Adjust concentrations for value boundary conditions

if ~isnan(A.sval)
    newtracer(1,:)=A.sval;
end

if ~isnan(A.bval)
    newtracer(max(Grd.nz),:)=A.bval;
end

%% Function VERTICALMIXING
function V = verticalmixing(tracer, Grd, mld, ent, In)
% VERTICALMIXING calculates vertical mixing of a tracer in the
% WCVIE2E_physicalmodel
%
% V = verticalmixing(tracer, Grd, mld, ent, In)
%
% This function is derived from Ianson and Allen. 2002. A two-dimensional
% nitrogen and carbon flux model in a coastal upwelling region. 
% Vertical mixing (tracer.d^-1) is estimated from the mixed layer depth, entrainment
% rates and other physical parameters.
%
%   Input variables:
%
%       tracer:  nz x nx array. Tracer concentrations in each spatial box at the
%                current time step (units vary based on tracer)
%
%       Grd   :  Structure holding spatial and temporal grid data for
%                WCVIE2E_physicalmodel simulations
%
%       mld   :  1 x nx array. Mixed layer depth at the current time step (m).
%
%       ent   :  1 x nx array. Entrainment rates at the current time step
%                (day^-1)
%
%       In    :  structure holding user-supplied input variables
%
%   Output variable:
%
%       V   :  nz x nx array. Flux(es) of tracer in the vertical dimension
%              that feed each spatial box (tracer.d^-1).

%Physical parameters

Mv=0.2;   % vertical mixing (m.d^-1)

% Creation of the matrix V

V=zeros(Grd.nz, Grd.nx);

%shelf & upper layer (i.e. mixed layer)
V(1,1)=((Mv./mld(1))+max(ent(1),0)).*(tracer(2,1)-tracer(1,1));

%shelf & lower layer
V(2,1)=((Mv./(Grd.zp(4,1)-In.dz(3,1)-mld(1)))+min(ent(1),0)).*(tracer(1,1)-tracer(2,1))...
    +(Mv./(Grd.zp(4,1)-In.dz(3,1)-mld(1))).*(tracer(3,1)-tracer(2,1));

%shelf & demersal layer
V(3,1)=(Mv./In.dz(3,1)).*(tracer(2,1)-tracer(3,1));

%slope & upper layer (i.e. mixed layer)
V(1,2)=((Mv./mld(2))+max(ent(2),0)).*(tracer(2,2)-tracer(1,2));

%slope & lower layer
V(2,2)=((Mv./(Grd.zp(4,2)-In.dz(3,2)-mld(2)))+min(ent(2),0)).*(tracer(1,2)-tracer(2,2))...
    +(Mv./(Grd.zp(4,2)-In.dz(3,2)-mld(2))).*(tracer(3,2)-tracer(2,2));

%slope & demersal layer
V(3,2)=(Mv./In.dz(3,2)).*(tracer(2,2)-tracer(3,2));


%% Function HORIZONTALMIXING
function H = horizontalmixing(tracer, Grd, mld, lb, In)
% HORIZONTALMIXING calculates horizontal mixing of a tracer in the
% WCVIE2E_physicalmodel
%
%H = horizontalmixing(tracer, Grd, mld, lateralforcing, In)
%
%This function is derived from Ianson and Allen. 2002. A two-dimensional
%nitrogen and carbon flux model in a coastal upwelling region. 
%
%
%   Input variables:
%
%       tracer:  nz x nx array. Tracer concentrations in each spatial box at the
%                current time step (units vary based on tracer)
%
%       Grd   : Structure holding spatial and temporal grid data for
%               WCVIE2E_physicalmodel simulations
%
%       mld   : 1 x nx array. Mixed layer depth at the current time step (m).
%
%       lb    : 1 x 5 array. Lateral boundary concentrations of tracer 
%               (units vary based on tracer). Column 1 = open ocean upper
%               layer; Column 2 = open ocean lower layer; Column 3 = rain;
%               Column 4 = freshwater from run-offs; Column 5 = VICC
%
%       In    :  structure holding user-supplied input variables
%
%
%   Output variable:
%
%       H   :  nz x nx array. Flux(es) of tracer in the horizontal dimension
%              that feed each spatial box (tracer.d^-1).

%Physical parameters

Mh=20;% horizontal mixing 

% Creation of the matrix H

H=zeros(Grd.nz, Grd.nx);

%shelf & upper layer (i.e. mixed layer)
H(1,1)=(Mh./In.dx(1,1)).*(tracer(1,2)-tracer(1,1));

%shelf & lower layer
H(2,1)=(Mh./In.dx(2,1)).*(tracer(2,2)-tracer(2,1));

%shelf & demersal layer
H(3,1)=(Mh./In.dx(3,1)).*(tracer(2,2)-tracer(3,1));

%slope & upper layer (i.e. mixed layer)
H(1,2)=((Mh./In.dx(1,2)).*(mld(1)./mld(2)).*(tracer(1,1)-tracer(1,2)))...
    +((Mh./In.dx(1,2)).*(lb(1)-tracer(1,2)));

%slope & lower layer
H(2,2)=((Mh./In.dx(2,2)).*((Grd.zp(4,1)-In.dz(3,1)-mld(1))./(Grd.zp(4,2)-In.dz(3,2)-mld(2))).*(tracer(2,1)-tracer(2,2)))...
    +((Mh./In.dx(2,2)).*(In.dz(3,1)./(Grd.zp(4,2)-In.dz(3,2)-mld(2))).*(tracer(3,1)-tracer(2,2)))...
    +((Mh./In.dx(2,2)).*(lb(2)-tracer(2,2)));

%slope & demersal layer
H(3,2)=(Mh./In.dx(3,2)).*(lb(2)-tracer(3,2));


%% Function ADVECTION

function X = advection(tracer, lb, In, Grd, mld, wndl, wndr, Sig, it)
% ADVECTION calculates advection of a tracer in the
% WCVIE2E_physicalmodel
%
%   X = advection(tracer, lb, In, Grd, mld, wndl, wndr, Sig, it)
%
%   This function is derived from Ianson and Allen. 2002. A two-dimensional
%   nitrogen and carbon flux model in a coastal upwelling region. 
%
%   Input variables:
%
%      tracer:  nz x nx array. Tracer concentrations in each spatial box at the
%               current time step (units vary based on tracer)
%
%      lb    :  1 x 5 array. Lateral boundary concentrations of tracer at the
%               current time step (units vary based on tracer). Column 1 =
%               upper open ocean; Column 2 = lower open ocean; Column 3 = rain;
%               Column 4 = freshwater from run-offs; Column 5 = VICC
%
%       Grd   :  Structure holding spatial and temporal grid data for
%                WCVIE2E_physicalmodel simulations
%
%       mld   :  1 x nx array. Mixed layer depth at the current time step (m).
%
%       In    :  structure holding user-supplied input variables
%
%       wndl  :  scalar. Local longshore wind stress (N.m^-2) at the
%                current time step.
%
%       wndr  : stucture related to remote wind forcing data
%
%       Sig   : nz x nx array holding the density (kg.m^-3) values in each
%               spatial box at the current time step
%
%       it    : current time step
%
%   Output variable:
%
%       X     :  nz x nx array. Flux(es) of tracer due to advection
%                that feed each spatial box (tracer.d^-1).
%
%
%
%-----------------------------
%Calculating Ekman velocity
%-----------------------------
%Physical parameters

cf= 2.*7.2921.*10.^(-5).*86400.*sin(In.Lat); %Coriolis parameter (rad.day.^-1)
Rd= 20000; %internal Rossby radius of deformation (m)

% calculating local a/d

loc=wndl./(mean2(Sig).*cf.*Rd);

if loc > 0
    loca=0;
    locd=loc; % northward wind (+ ve) causes downwelling (-ve)
else
    loca=-loc;
    locd=0; % southward wind (-ve) causes upwelling (+ve)
end

%--------------------------------------
%Calculating remote upwelling velocity
%--------------------------------------

%Evaluate coefficients of polynomial fit

tfil=wndr.t;
xfil=wndr.data;% wind stress based remote upwelling index

% When current time steps not equal to 1 or Grd.nt
    
%if it~= 1 || Grd.nt

%det = tfil(it-1)*tfil(it-1)*(tfil(it)-tfil(it+1))-tfil(it)*tfil(it)* ...
%      (tfil(it-1)-tfil(it+1)) + tfil(it+1)*tfil(it+1)*(tfil(it-1)-tfil(it));

%aa = ((tfil(it)-tfil(it+1))*xfil(it-1)-(tfil(it-1)-tfil(it+1))*xfil(it)+ ...
%     (tfil(it-1)-tfil(it))*xfil(it+1))./det;

%bb = (-1*(tfil(it)*tfil(it)-tfil(it+1)*tfil(it+1))*xfil(it-1)+(tfil(it-1)*tfil(it-1)- ...
%     tfil(it+1)*tfil(it+1))*xfil(it)-(tfil(it-1)*tfil(it-1)-tfil(it)*tfil(it))*xfil(it+1))./det;

%cc = (xfil(it-1)*(tfil(it)-tfil(it+1))*tfil(it+1)*tfil(it) ...
%     -xfil(it)*(tfil(it-1)-tfil(it+1))*tfil(it+1)*tfil(it-1) ...
%     +xfil(it+1)*(tfil(it-1)-tfil(it))*tfil(it)*tfil(it-1))./det;

% when current time step = 1
%elseif it==1
    
%det = (-1)*(-1)*(tfil(it)-tfil(it+1))-tfil(it)*tfil(it)*(-1-tfil(it+1))+...
%      tfil(it+1)*tfil(it+1)*(-1-tfil(it));

%aa = ((tfil(it)-tfil(it+1))*0.951724-(-1-tfil(it+1))*xfil(it)+ ...
%     (0-tfil(it))*xfil(it+1))./det;
 
%bb = (-1*(tfil(it)*tfil(it)-tfil(it+1)*tfil(it+1))*0.951724+(-1*(-1)- ...
%     tfil(it+1)*tfil(it+1))*xfil(it)-(-1*(-1)-tfil(it)*tfil(it))*xfil(it+1))./det;
 
%cc = (0.951724*(tfil(it)-tfil(it+1))*tfil(it+1)*tfil(it) ...
%      -xfil(it)*(-1-tfil(it+1))*tfil(it+1)*(-1) ...
%      +xfil(it+1)*((-1)-tfil(it))*tfil(it)*(-1))./det;
 
% When current time step = Grd.nt
%else
    
%det= tfil(it-1)*tfil(it-1)*(tfil(it)-(Grd.nt+1))-tfil(it)*tfil(it)* ...
%      (tfil(it-1)-(Grd.nt+1)) + (Grd.nt+1)*(Grd.nt+1)*(tfil(it-1)-tfil(it));

%aa = ((tfil(it)-(Grd.nt+1))*xfil(it-1)-(tfil(it-1)-(Grd.nt+1))*xfil(it)+ ...
%     (tfil(it-1)-tfil(it))*1.35953)./det;

%bb = (-1*(tfil(it)*tfil(it)-(Grd.nt+1)*(Grd.nt+1))*xfil(it-1)+(tfil(it-1)*tfil(it-1)- ...
%     (Grd.nt+1)*(Grd.nt+1))* xfil(it)-(tfil(it-1)*tfil(it-1)-tfil(it)*tfil(it))*1.35953)./det;

%cc = (xfil(it-1)*(tfil(it)-(Grd.nt+1))*(Grd.nt+1)*tfil(it) ...
%     -xfil(it)*(tfil(it-1)-(Grd.nt+1))*(Grd.nt+1)*tfil(it-1) ...
%     +1.35953*(tfil(it-1)-tfil(it))*tfil(it)*tfil(it-1))./det;
% 
%end
    
%evaluate force 
%force=aa*Grd.time(it)*Grd.time(it)+bb*Grd.time(it)+cc;

if xfil(it) > 0
    rema=0;
    remd=xfil(it);%northward wind (+ve) causes downwelling (-ve)
else
    rema=-xfil(it);
    remd=0;%southward wind (-ve) causes upwelling (+ve)
end

% add together loca/d and rema/d

a=loca+rema;
d=locd+remd;

afl=a.*In.dx(1,2);%a*slope width
dfl=d.*In.dx(1,2);%d*slope width

%creating X vector
X=zeros(Grd.nz, Grd.nx);

%%Buoyancy fluxes

%flux per length of coastline for rainfall

P=0.03.*exp(-0.3.*(5-2.*cos(2.*pi.*(Grd.time(it)+20)/365)))...
   -0.03.*exp(-0.6.*(10-6.*cos(2.*pi/365.*(Grd.time(it)-150))));

%flux per length of coastline for terrigenous runoff

R=0.5.*(0.11.*exp(-0.3.*(5-2.*cos(2.*pi.*(Grd.time(it)+20)/365)))...
   -0.11.*exp(-0.6.*(10-6.*cos(2.*pi/365.*(Grd.time(it)-150)))));

%flux per length of coastline for the VICC

C=0.5.*exp(-0.6.*(5-2.*cos(2.*pi.*(Grd.time(it)-150)/365)));

%Buoyancy fluxes

Bsh=((P./mld(1)).*(lb(3)-tracer(1,1)))+((R./mld(1)).*(lb(4)-tracer(1,1)))...
    +((C./mld(1)).*(lb(5)-tracer(1,1)));%shelf box buoyancy fluxes

Bsl=((P./mld(1)).*(lb(3)-tracer(1,1)));%slope box buoyancy fluxes


%shelf & upper layer (i.e. mixed layer)

X(1,1)= (((afl./In.dx(1,1))./mld(1)).*(tracer(2,1)-tracer(1,1)))...
    +(((dfl./In.dx(1,1))./mld(1)).*(tracer(1,2)-tracer(1,1)))+Bsh;

%shelf & lower layer

X(2,1)=(((afl./In.dx(2,1))./(Grd.zp(4,1)-In.dz(3,1)-mld(1))).*(tracer(2,2)-tracer(2,1)))...
    +(((dfl./In.dx(2,1))./(Grd.zp(4,1)-In.dz(3,1)-mld(1))).*(tracer(1,1)-tracer(2,1)));

%shelf & demersal layer: no upwelling advection occuring in this box

X(3,1)=0;

%slope & upper layer (i.e. mixed layer)

X(1,2)=(((afl./In.dx(1,2))./mld(2)).*(tracer(1,1)-tracer(1,2)))...
    +(((dfl./In.dx(1,2))./mld(2)).*(lb(1)-tracer(1,2)))+Bsl;

%slope & lower layer

X(2,2)=(((afl./In.dx(2,2))./(Grd.zp(4,2)-In.dz(3,2)-mld(2))).*(lb(2)-tracer(2,2)))...
    +(((dfl./In.dx(2,2))./(Grd.zp(4,2)-In.dz(3,2)-mld(2))).*(tracer(2,1)-tracer(2,2)));

%slope & demersal layer
    
X(3,2)=0;
