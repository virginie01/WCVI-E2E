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









