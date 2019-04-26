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