function [dxdt, V, H, X, varargout] = mixtracerTS_ode(time, tracer, In, Grd, mld, ...
    ent, lb, Sig, tauy, prate, varargin)
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
%      mld     :  Structure holding mixed layer depth data (m). Data:
%                 datant+1 x nx array.
%
%      ent     :  structure holding entrainment rate data (s^-1). Data:
%                 datant+1 x nx array.
%
%      lb      :  Structure holding lateral boundary concentrations of 
%                 tracer (units vary based on tracer).
%                 Data: datant+1 x 6 array. Column 1 = open ocean upper layer;
%                 Column 2 = open ocean lower layer; Column 3 = rain shelf;
%                 Column 4: rain slope; Column 5 = freshwater from run-offs; 
%                 Column 6 = VICC
%
%      Sig     :  nz x nx array holding the density (kg.m^-3) values in each
%                 spatial box at the current time step
%
%      tauy    :  Structure holding N-S wind stress data (N.m^-2).
%                 Data: datant+1 x nx array.
%                 datant+1 x 1 array if from tauy2
%
%      prate   :  Structure holding precipitation rates. Expressed 
%                 in m.s^-1. Data: datant+1 x nx array 
%
%
%   Optional input variables passed as parameter + structures needed to
%   calculate the parameter (i.e. relevant for temperature):
%
%      sflux   :   heat fluxes across the surface. 1 x nx array. Expressed 
%                  in W.m-2                 
%
%      bflux   :   heat fluxes across the seafloor. 1 x nx array. Expressed
%                  in W.m-2
%
%      source  :   heat fluxes throughout the water column in all boxes. 
%                  nz x nx array. Expressed in W.m-2
%
%    Output variable:
%
%      dxdt    :   nz x nx array. ODE for tracer concentrations(tracer unit.s^-1).

%% Calculating each component of mixing and advection

V = verticalmixing(tracer, Grd, mld, ent, In, time);


H = horizontalmixing(tracer, Grd, mld, lb, In, time);


X = advection(tracer, lb, In, Grd, mld, tauy, Sig, prate, time);


%% Check input

if ~isequal(size(tracer), size(V), size(H), size (X))
    error('flux coefficient array must be same size as tracer array');
end

%% Calculating new tracer concentrations
dxdt = V+H+X;

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
        [~, qs, ql, qlw] = WCVIE2E_calcheat(qi, airtmp, dewptT, tracer(1,:),...
            Wspd10, Qo, Grd.nx);
    % Sensible, latent, and longwave fluxes are applied at the surface
    
    srfhflx = qs + ql + qlw;
    
    % Heat fluxes (W/m^2) are converted to temperature fluxes (deg C/s)
    
    Cp = 4125; % Typical specific heat for water (J kg^-1 degC^-1) @10 deg., 32 ppt)
    srf_Tflx = srfhflx./(Sig(1,:).*Cp.* mlayer);
    
    % first optional output 
    
    varargout{1}=srf_Tflx;
    
    % calculate new ode including surface heat fluxes
    dxdt(1,:)= dxdt(1,:)+ srf_Tflx;
    
    end
    
    % source heat fluxes
    
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
       
        zp = [Grd.zp(1,1), Grd.zp(1,2); ...
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
        solhflx(:,i) = (1-In.alb) .* qi(i) * diff(...
        In.prad1*exp(In.krad1*zp(:,i)) + (1-In.prad1)*exp(In.krad2*zp(:,i)));        
        end
        
    % Heat fluxes (W/m^2) are converted to temperature fluxes (deg C/s)
    
    Cp = 4125; % Typical specific heat for water (J kg^-1 degC^-1) @10 deg., 32 ppt)
    
    sol_Tflx = solhflx./(Sig.* Cp.* dz);
    
    % second optional output
    
    varargout{2} = sol_Tflx;
    
    % new tracer concentration
    dxdt = dxdt + sol_Tflx;

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
%       mld   :  Structure holding mixed layer depth data (m). Data:
%                datant+1 x nx array.
%
%       ent   :  Structure holding entrainment rate data (s^-1). Data:
%                datant+1 x nx array
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

Mv=0.2./86400;   % vertical mixing (m.s^-1)

% Creation of the matrix V

V=zeros(Grd.nz, Grd.nx);

% finding time index corrsponding to current time

idx=mld.t==time & ent.t == time;

%shelf & upper layer (i.e. mixed layer)
V(1,1)=((Mv./mld.data(idx,1))+max(ent.data(idx,1),0)).*(tracer(2,1)-tracer(1,1));

%shelf & lower layer
V(2,1)=((Mv./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))+min(ent.data(idx,1),0)).*(tracer(1,1)-tracer(2,1))...
    +(Mv./(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1))).*(tracer(3,1)-tracer(2,1));

%shelf & demersal layer
V(3,1)=(Mv./In.dz(3,1)).*(tracer(2,1)-tracer(3,1));

%slope & upper layer (i.e. mixed layer)
V(1,2)=((Mv./mld.data(idx,2))+max(ent.data(idx,2),0)).*(tracer(2,2)-tracer(1,2));

%slope & lower layer
V(2,2)=((Mv./(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2)))...
    + min(ent.data(idx,2),0)).*(tracer(1,2)-tracer(2,2))...
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
%       mld   : Structure holding mixed layer depth data (m). Data:
%               datant+1 x nx array.
%
%       lb    : Structure holding lateral boundary concentrations of tracer 
%               (units vary based on tracer). Data: datant+1 x 6 array. 
%               Column 1 = open ocean upper layer; 
%               Column 2 = open ocean lower layer; Column 3 = rain shelf;
%               Column 4: rain slope; Column 5 = freshwater from run-offs; 
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

cf = 2*7.2921*10.^(-5)*sin(In.Lat); %Coriolis parameter (rad.s^-1)
Rd = 20000; %internal Rossby radius of deformation (m)

% calculating local a/d

uek=tauy.data(idx)/(Sig(1,1)*cf*Rd); %uek in m.s^-1

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

%P= 0.03.*exp(-0.3.*(5-2.*cos(2.*pi.*(Grd.time(it)+20)/365)))...
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
% Bsh=0;

%slope box buoyancy fluxes
Bsl=((prate.data(idx,2)./mld.data(idx,2)).*(lb.data(idx,4)-tracer(1,2)));
%Bsl=0;

%shelf & upper layer (i.e. mixed layer)

X(1,1)= (((afl/In.dx(1,1))/mld.data(idx,1))*(tracer(2,1)-tracer(1,1)))...
    +(((dfl/In.dx(1,1))/mld.data(idx,1))*(tracer(1,2)-tracer(1,1)))+Bsh;

%shelf & lower layer

X(2,1)=(((afl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(2,2)-tracer(2,1)))...
    +(((dfl/In.dx(2,1))/(Grd.zp(4,1)-In.dz(3,1)-mld.data(idx,1)))*(tracer(1,1)-tracer(2,1)));

%shelf & demersal layer: no upwelling advection occuring in this box

X(3,1)=0;

%slope & upper layer (i.e. mixed layer)

X(1,2)=(((afl/In.dx(1,2))/mld.data(idx,2))*(tracer(1,1)-tracer(1,2)))...
    +(((dfl/In.dx(1,2))/mld.data(idx,2))*(lb.data(idx,1)-tracer(1,2)))+Bsl;

%slope & lower layer

X(2,2)=(((afl/In.dx(2,2))/(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2)))*(lb.data(idx,2)-tracer(2,2)))...
    +(((dfl/In.dx(2,2))/(Grd.zp(4,2)-In.dz(3,2)-mld.data(idx,2)))*(tracer(2,1)-tracer(2,2)));

%slope & demersal layer
    
X(3,2)=0;
