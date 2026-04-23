function [newtracer,db,Flx] = verticalflux(tracer, wsink, dt, G)
%VERTICALFLUX Calculates vertical movement within mixed layer model
%
% newtracer = verticalflux(tracer, wsink, dt, dz, openbot)
%
% Calculated changes in tracer concentration due to non-mixing processes.
%
% Input variables:
%
%   tracer:     nz x nx x nbsv array, concentration of tracer
%
%   wsink:      nz x nx x nbsv array, vertical velocities (m/s)
%
%   dt:         scalar, model time step (s)
%
%   G:          Structure containing grid parameters for current time step
%
% Output variables:
%
%   newtracer:  nz x nx x nbsv array, new value of tracer concentrations vertical flux over the 
%               current time interval

% Copyright 2009 Kelly Kearney

% Calculate amount of tracer leaving each depth layer

if any(tracer < 0)
    error('Negative tracer');
end

[~, ~, nb] = size(tracer);

V = G.dz.*G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2)));% nz x nx m3
V = repmat(V,1,1,nb);%nz x nx x nb m3

dz = repmat(G.dz,1,1,nb);% nz x nx x nb

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

% No flux through ocean surface or floor (unless open bottom)
 a = fluxout(1,:,:)>0;
 fluxout(a)=0;
 b = fluxout(end,:,:)<0;
 fluxout(b)=0;

% Characterize flux out of each layer as going up or down

fluxdown = zeros(size(tracer));
fluxup   = zeros(size(tracer));

isup = fluxout > 0;
fluxup(isup) = fluxout(isup);
fluxdown(~isup) = -fluxout(~isup);% flux down were initially negative and now positive

% Calculate total change due to loss from a layer and gain from adjacent
% layers

tracerchange = zeros(size(tracer)); % nz x nx x nbsv mol N
tracerchange(2:end-1,:,:) = fluxup(3:end,:,:) + fluxdown(1:end-2,:,:) - abs(fluxout(2:end-1,:,:));
tracerchange(1,:,:)       = fluxup(2,:,:)                             - abs(fluxout(1,:,:));
tracerchange(end,:,:)     =                     fluxdown(end-1,:,:)   - abs(fluxout(end,:,:));

tracerchange = tracerchange./V; % nz x nx x nb molN/m3 leaving each box

db = tracerchange./dt; %nz x nx x nbsv molN.m-3.s-1

[nz, nx, nbsv] = size(tracer);
Flx.vflx = zeros(nbsv+2, nbsv+2, nz, nx);

for i = 1:nbsv
Flx.vflx(i,i,:,:) = permute(db(:,:,i), [3 4 1 2]);
end

newtracer = tracer + tracerchange;
%newtracer = tracer;
