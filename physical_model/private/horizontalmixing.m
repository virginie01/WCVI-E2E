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
