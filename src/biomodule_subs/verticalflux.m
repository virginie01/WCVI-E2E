function [newtracer,db,Flx] = verticalflux(tracer, wsink, dt, G)
% VERTICALFLUX Simulates vertical sinking/rising of biological tracers.
%
%   [newtracer, db, Flx] = verticalflux(tracer, wsink, dt, G)
%
%   This function computes the vertical flux of tracers (e.g., plankton)
%   due to sinking or buoyancy-driven movement using a layer-based approach.
%
%   INPUTS:
%       tracer   : [nz x nx x nbsv] array of tracer concentrations (mol N/m³)
%       wsink    : [nz x nx x nbsv] array of vertical velocities (m/s)
%                  Positive = upward, Negative = downward
%       G        : Structure containing grid parameters for current time step
%
%   OUTPUTS:
%       newtracer: [nz x nx x nbsv] array of updated tracer concentration after vertical flux
%       db       : tracer rate of change [mol N/m³/s]
%       Flx      : structure with vertical flux diagnostics
%                   - vflx: [nbsv+2 x nbsv+2 x nz x nx] flux matrix
%
%   Notes:
%       - No flux through ocean surface or bottom layer.
%       - Assumes closed top/bottom boundary unless otherwise handled externally.
%       - Tracers are conserved vertically (redistributed, not lost).
%
% This file was initially derived from the original verticalflux
% routine developed by Kelly Kearney for the WCE/NEMURO framework
% and substantially extended for the WCVI-E2E coastal upwelling
% ecosystem model.
%
% Original framework:
% Copyright (c) 2009 Kelly Kearney
%
% Major modifications and extensions by Virginie Bornarel (2017–2026)
% include:
%   - adaptation from 1D to 2D multi-tracer transport
%   - explicit grid-cell volume scaling using WCVI geometry
%   - support for nz x nx x nbsv biological state arrays
%   - implementation of full-box transfer handling for large vertical velocities
%   - expanded vertical flux diagnostics and bookkeeping
%   - integration with WCVI-E2E transport and ecosystem flux structures
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% --------- INITIAL SAFETY CHECKS ---------
if any(tracer < 0)
    error('Negative tracer');
end

[~, ~, nb] = size(tracer);

% --------- VOLUME CALCULATION (m³ per grid cell) ---------
V = G.dz.*G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2)));% nz x nx m3
V = repmat(V,1,1,nb);%nz x nx x nb m3

dz = repmat(G.dz,1,1,nb);% nz x nx x nb

% --------- INITIALIZE FLUXES ---------
fluxout = zeros(size(tracer));

% if the distance travelled during a time step is less than dz, just a
% proportion of ZL2 leaves the box
p = abs(wsink.*(dt./dz))<1; 
fluxout(p) = tracer(p).*wsink(p).*(dt./dz(p));% nz x nx x nb, molN/m3
fluxout(p) = fluxout(p).*V(p);% nz x nx x nb, mol N leaving each box

%if distance travelled within a time step is bigger than dz then all ZL2
%leave the box
allpos = (wsink.*(dt./dz))>=1;
fluxout(allpos) = tracer(allpos).*V(allpos);%nz x nx x nb mol N leaving each box, + for upward movement

allneg = (wsink.*(dt./dz))<=-1;
fluxout(allneg) = -tracer(allneg).*V(allneg);%nz x nx x nb mol N leaving each box, - for downward movement

% --------- BLOCK TOP AND BOTTOM FLUXES ---------
 a = fluxout(1,:,:)>0;
 fluxout(a)=0;
 b = fluxout(end,:,:)<0;
 fluxout(b)=0;


% --------- DECOMPOSE INTO UPWARD AND DOWNWARD FLUX ---------
fluxdown = zeros(size(tracer));
fluxup   = zeros(size(tracer));

isup = fluxout > 0;
fluxup(isup) = fluxout(isup);
fluxdown(~isup) = -fluxout(~isup);% flux down were initially negative and now positive

% --------- CALCULATE TRACER CHANGES ---------
tracerchange = zeros(size(tracer)); % nz x nx x nbsv mol N
tracerchange(2:end-1,:,:) = fluxup(3:end,:,:) + fluxdown(1:end-2,:,:) - abs(fluxout(2:end-1,:,:));
tracerchange(1,:,:)       = fluxup(2,:,:)                             - abs(fluxout(1,:,:));
tracerchange(end,:,:)     =                     fluxdown(end-1,:,:)   - abs(fluxout(end,:,:));

tracerchange = tracerchange./V; % nz x nx x nb molN/m3 leaving each box

db = tracerchange./dt; %nz x nx x nbsv molN.m-3.s-1

% --------- CREATE FLUX STRUCTURE ---------
[nz, nx, nbsv] = size(tracer);
Flx.vflx = zeros(nbsv+2, nbsv+2, nz, nx);

for i = 1:nbsv
    Flx.vflx(i,i,:,:) = permute(db(:,:,i), [3 4 1 2]);
end

% --------- UPDATE TRACER FIELD ---------
newtracer = tracer + tracerchange;
%newtracer = tracer;
