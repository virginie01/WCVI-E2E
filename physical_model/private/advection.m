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
%      Grd   :  Structure holding spatial and temporal grid data for
%               WCVIE2E_physicalmodel simulations
%
%      mld   :  1 x nx array. Mixed layer depth at the current time step (m).
%
%      In    :  structure holding user-supplied input variables
%
%      wndl  :  scalar. Local longshore wind stress (N.m^-2) at the
%               current time step.
%
%      wndr  : stucture related to remote wind forcing data
%
%      Sig   : nz x nx array holding the density (kg.m^-3) values in each
%              spatial box at the current time step
%
%      it    : current time step
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
    loca=0
    locd=loc % northward wind (+ ve) causes downwelling (-ve)
else
    loca=-loc
    locd=0 % southward wind (-ve) causes upwelling (+ve)
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

if xfil > 0
    rema=0;
    remd=xfil;%northward wind (+ve) causes downwelling (-ve)
else
    rema=-xfil;
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










