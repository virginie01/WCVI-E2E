function X = advection(tracer, lb, In, Grd, mld, tauy, Sig, prate, time)
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
%      lb    :  Structure holding lateral boundary concentrations of tracer 
%               (units vary based on tracer). Data: datant+1 x 6 array. 
%               Column 1 = open ocean upper layer; 
%               Column 2 = open ocean lower layer; Column 3 = rain shelf;
%               Column 4: rain slope; Column 5 = freshwater from run-offs; 
%               Column 6 = VICC
%
%      In    :  Structure holding user-supplied input variables
%
%      Grd   :  Structure holding spatial and temporal grid data for
%               WCVIE2E_physicalmodel simulations
%
%      mld   :  Structure holding mixed layer depth data (m). Data:
%               datant+1 x nx array.
%
%      tauy  :  Structure holding N-S wind stress data (N.m^-2).
%               Data: datant+1 x nx array.
%               datant+1 x 1 array if from tauy2
%
%      Sig   :  nz x nx array holding the density (kg.m^-3) values in each
%               spatial box at the current time-step
%
%      prate :  Structure holding precipitation rates. Expressed 
%               in m.s^-1. Data: datant+1 x nx array 
%
%      time  :  Scalar. Time elapsed since start of the simulation (s).
%
%   Output variable:
%
%       X     :  nz x nx array. Flux(es) of tracer due to advection
%                that feed each spatial box (tracer.s^-1) at the current
%                time step
%
% logical vector indicating current time 

 idx = lb.t==time & mld.t==time & tauy.t==time & prate.t==time;

%-----------------------------
%Calculating Ekman velocity
%-----------------------------
%Physical parameters

cf= 2.*7.2921.*10.^(-5).*sin(In.Lat); %Coriolis parameter (rad.s^-1)
Rd= 20000; %internal Rossby radius of deformation (m)

% calculating local a/d

uek=tauy.data(idx,1)./(Sig(1,1).*cf.*Rd); %uek in m.s^-1

if uek > 0
    a=0;
    d=uek; % northward wind (+ ve) causes downwelling (-ve)
else
    a=-uek;
    d=0; % southward wind (-ve) causes upwelling (+ve)
end

%--------------------------------------
%Calculating remote upwelling velocity
%--------------------------------------

%Evaluate coefficients of polynomial fit

%tfil=wndr.t;
%xfil=wndr.data;% wind stress based remote upwelling index

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

%if xfil(it) > 0
%    rema=0;
%    remd=xfil(it);%northward wind (+ve) causes downwelling (-ve)
%else
%    rema=-xfil(it);
%    remd=0;%southward wind (-ve) causes upwelling (+ve)
%end

% add together loca/d and rema/d

%a=loca+rema;
%d=locd+remd;

afl=a.*In.dx(1,2);%a*slope width
dfl=d.*In.dx(1,2);%d*slope width

%creating X vector
X=zeros(Grd.nz, Grd.nx);

%%Buoyancy fluxes

%flux per length of coastline for rainfall

%P=0.03.*exp(-0.3.*(5-2.*cos(2.*pi.*(Grd.time(it)+20)/365)))...
%   -0.03.*exp(-0.6.*(10-6.*cos(2.*pi/365.*(Grd.time(it)-150))));

%flux per length of coastline for terrigenous runoff: original cycle from
%Debby's publication in m.day-1.

R=0.05.*exp(-0.5.*(3.5-(2.*cos((2.*pi./365).*((time./86400)+20)))));
R=R./86400;

%flux per length of coastline for the VICC: original cycle from
%Debby's publication in m.day-1.

C=0.05.*exp(-0.6.*(5-(2.*cos((2.*pi./365).*((time./86400)+150)))));
C=C./86400;

%Buoyancy fluxes

%shelf box buoyancy fluxes
Bsh=((prate.data(idx,1)./mld.data(idx,1)).*(lb.data(idx,3)-tracer(1,1)))...
    +((R./mld.data(idx,1)).*(lb.data(idx,5)-tracer(1,1)))...
    +((C./mld.data(idx,1)).*(lb.data(idx,6)-tracer(1,1)));

%slope box buoyancy fluxes
Bsl=((prate.data(idx,2)./mld.data(idx,2)).*(lb.data(idx,4)-tracer(1,2)));


%shelf & upper layer (i.e. mixed layer)

X(1,1)= (((afl./In.dx(1,1))./mld.data(idx,1)).*(tracer(2,1)-tracer(1,1)))...
    +(((dfl./In.dx(1,1))./mld.data(idx,1)).*(tracer(1,2)-tracer(1,1)))+Bsh;

%shelf & lower layer

X(2,1)=(((afl./In.dx(2,1))./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(tracer(2,2)-tracer(2,1)))...
    +(((dfl./In.dx(2,1))./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(tracer(1,1)-tracer(2,1)));

%shelf & demersal layer: no upwelling advection occuring in this box

X(3,1)=0;

%slope & upper layer (i.e. mixed layer)

X(1,2)=(((afl./In.dx(1,2))./mld.data(idx,2)).*(tracer(1,1)-tracer(1,2)))...
    +(((dfl./In.dx(1,2))./mld.data(idx,2)).*(lb.data(idx,1)-tracer(1,2)))+Bsl;

%slope & lower layer

X(2,2)=(((afl./In.dx(2,2))./(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2))).*(lb.data(idx,2)-tracer(2,2)))...
    +(((dfl./In.dx(2,2))./(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2))).*(tracer(2,1)-tracer(2,2)));

%slope & demersal layer
    
X(3,2)=0;
