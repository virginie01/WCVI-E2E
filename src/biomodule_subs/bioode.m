function [db, Flx, Diag] = bioode(nemflag,~, bio, G, P, B, Arch, it)
%BIOODE Main ODE RHS for WCVI E2E module (or standalone NEMURO model if
%nemflag is turned on)
%
% [db, Flx, Diag] = bioode(nemflag, time, bio, G, P, B, Arch, it)
%
% PURPOSE
%   Computes source/sink terms and derivatives for all biological state 
%   variables in each (z,x) box, suitable for MATLAB ODE solvers.
%
% INPUTS
%   nemflag : logical scalar
%       true  -> NEMURO standalone version
%       false -> full WCVI-E2E model
%
%   bio : [nz nx nb] array
%       state variable concentrations at beginning of time step
%
%   G : struct
%       Grid geometry
%
%   P : struct
%       Physical forcing
%
%   B : struct
%       Biological parameters, indexing (B.idx.*), trophic link codes (B.links),
%       import fractions (B.import), nekton mask (B.isnek), diapause flags, etc.
%
%   Arch : Struct array for archiving diagnostic outputs
% 
%   it : Current simulation time index
%
% OUTPUTS
%   db : [nz x nx x nb x 1] numeric
%       dB/dt (mol N m^-3 s^-1) for each state variable in each box.
%
%   Flx : struct
%       Fluxes between sources/sinks in each box. Fields are 4-D arrays 
%       [nb+2 x nb+2 x nz x nx]
%
%   Diag : struct
%       Intermediate fluxes for diagnostic purposes

bio = max(bio, 0); % Negative biomass treated as 0

%------------------------------
% Set up various bio arrays
%------------------------------
if B.diapause    
   [bv, ~ , ~, basum1, basum2, bfrac, nb, nz, nx, zlfrac] =  biomasssetup(bio, G, B);
else
   [bv, ~ , ~, basum1, basum2, bfrac, nb, nz, nx] =  biomasssetup(bio, G, B);   
end

%------------------------------
% Photosynthesis-related
% fluxes (gpp, exc, resp)
%------------------------------

[gpp, exc, resp, ~, psmax, Lfc, no3lim, nh4lim, silim, I, ...
    kappa, kappaP] = primprod(bv.orig, G, P, B, nz, nx, nb);


%------------------------------
% Aydin Ecosim primary 
% production
%------------------------------
if ~nemflag
if B.ecosimppflag  
    npp = zeros(nb+2,nb+2,nz,nx);
    istoodeep = -G.z < -100;
    for ib = [B.idx.ps B.idx.pl]
        h = exp(-2) + 1;  % TODO
        ybio = bv.orig(:,:,ib) ./ B.b0v(ib);
        npp(B.idx.mys,ib,:,:) = B.p0v(ib) .* (h.*ybio)./(h-1+ybio);
        npp(B.idx.mys,ib,istoodeep,:) = 0;
    end
end
end

%------------------------------
% Predation (involves nekton)
%------------------------------

if ~nemflag 
    
pred_out = zeros(nb+2,nb+2,nz,nx);
pred_in = zeros(nb+2,nb+2,nz,nx);%this one is to accommodate the N-Z predation fluxes that enter nekton groups in demersal layer based on pred1
 
% Predation flux based on aydin ecosim functional response
% (ZL never predate, so don't have to worry about separating in/out here)
% pred1 for pred_in based on true densities/pred2 for pred_out based on
% concentrations

pred1 = aydinfrnew(basum1.prey, basum1.pred, B.b0, B.q0, B.x, B.d, B.theta); % nb x nb, mol N m^-2 s^-1
pred2 = aydinfrnew(basum2.prey, basum2.pred, B.b0, B.q0, B.x, B.d, B.theta); % nb x nb, mol N m^-2 s^-1

% For nekton-nekton links, this total flux sits in box 31

prednn = zeros(nb);
prednn(B.links == 3) = pred1(B.links == 3);
prednn = (prednn.*G.area)./(G.dz(3,1).*G.dx(3,1).*(G.area./(G.dx(1,1)+G.dx(1,2)))); % nb x nb, mol N m^-3 s^-1

% For nekton-zooplankton links and for the zooplankton side, the flux is distributed based on where the
% prey was located (based on pred2)- For the nekton side, the flux is located in demersal layer and based on pred1 

prednz1 = zeros(nb);
prednz1(B.links == 2 | B.links == 5 | B.links == 7) = pred1(B.links == 2 | B.links == 5 | B.links == 7);
prednz1 = (prednz1.*G.area)./(G.dz(3,1).*G.dx(3,1).*(G.area./(G.dx(1,1)+G.dx(1,2)))); % nb x nb, mol N m^-3 s^-1

prednz2 = pred2.*G.area; %nb x nb, molN/s
prednz2 = bsxfun(@times, prednz2, bfrac.prey); %nb x nb x nz x nx molN/s
vpb = zeros(nb,nb,nz,nx); %nb x nb x nz x nx
for i =1:nb
    for j =1:nb
vpb(i,j,:,:)=G.dz.*G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2)));
    end
end

prednz2 = prednz2./vpb; %nb x nb x nz x nx molN/m3/s
prednz2(isnan(prednz2)) = 0;

for iprey = 1:nb
    for ipred = 1:nb
        if B.links(iprey,ipred) ~= 2 && B.links(iprey,ipred) ~= 5 &&...
               B.links(iprey,ipred) ~= 7
           prednz2(iprey,ipred,:,:) = 0;
        end
    end
end

% Combine

pred_out(1:nb,1:nb,3,1) = prednz2(:,:,3,1) + prednn;
pred_out(1:nb,1:nb,1,1) = prednz2(:,:,1,1); pred_out(1:nb,1:nb,2,1) = prednz2(:,:,2,1);
pred_out(1:nb,1:nb,1,2) = prednz2(:,:,1,2);pred_out(1:nb,1:nb,2,2) = prednz2(:,:,2,2);
pred_out(1:nb,1:nb,3,2) = prednz2(:,:,3,2);

pred_in(1:nb,1:nb,3,1) = prednz1 + prednn;

% Include imports

predimp = bsxfun(@times, bv.orig, reshape(B.import,1,1,nb)); %nz x nx x nb, molN/m3/s

pred_out(B.idx.mys,1:nb,:,:) = permute(predimp, [4 3 1 2]);
pred_out(B.idx.mys,[~B.isnek' false false],:,:)=0;

pred_in(B.idx.mys,1:nb,3,1) = reshape(predimp(3,1,:),1,nb,1,1);
pred_in(B.idx.mys,[~B.isnek' false false],3,1)=0;

else
    
pred_out = zeros(nb+2,nb+2,nz,nx);


end

%------------------------------
% Grazing (involves only 
% plankton, and includes temp
% influence)
%------------------------------

if nemflag

q = zeros (nb+2,nb+2,nz,nx);

for iz = 1:nz
    for ix = 1:nx
    q(1:nb,1:nb,iz,ix) = ivlev(bv.prey(:,:,iz,ix), bv.pred(:,:,iz,ix),...
        B.grmax, B.Kgra, B.lambda, B.thresh, P.T(iz,ix)); % mol N m^-3 s^-1
    end    
end

q(B.idx.pl,B.idx.zp,:,:)= q(B.idx.pl,B.idx.zp,:,:).*reshape(exp(-B.grpusai(2,5).*(bv.zlcombo(:,:,B.idx.zl)+...
    bv.zlcombo(:,:,B.idx.zs))),1,1,nz,nx);

q(B.idx.zs,B.idx.zp,:,:)= q(B.idx.zs,B.idx.zp,:,:).*reshape(exp(-B.grpusai(3,5).*bv.zlcombo(:,:,B.idx.zl)),1,1,nz,nx);

graze = max(0,q);

elseif ~nemflag

q = zeros (nb+2,nb+2,nz,nx);
for iz = 1:nz
    for ix =1:nx
    q(1:nb,1:nb,iz,ix) = ivlev(bv.prey(:,:,iz,ix), bv.pred(:,:,iz,ix),...
        B.grmax, B.Kgra, B.lambda, B.thresh, P.T(iz,ix)); % mol N m^-3 s^-1
    end
end

q(B.idx.pl,B.idx.zp,:,:)= q(B.idx.pl,B.idx.zp,:,:).*reshape(exp(-B.grpusai(2,5).*(bv.zlcombo(:,:,B.idx.zl)+...
    bv.zlcombo(:,:,B.idx.zs))),1,1,nz,nx);

q(B.idx.zs,B.idx.zp,:,:)= q(B.idx.zs,B.idx.zp,:,:).*reshape(exp(-B.grpusai(3,5).*bv.zlcombo(:,:,B.idx.zl)),1,1,nz,nx);

graze = max(0,q);

% Involves only plankton

for iprey = 1:nb
    for ipred = 1:nb
        if B.links(iprey,ipred) ~= 1 && B.links(iprey,ipred) ~= 4 &&...
               B.links(iprey,ipred) ~= 6 
           graze(iprey,ipred,:,:) = 0;
        end
    end
end

% Include imports

grazeimp = bsxfun(@times, bv.orig, reshape(B.import,1,1,nb)); %nz x nx x nb, molN/m3/s
graze(B.idx.mys,1:nb,:,:) = permute(grazeimp, [4 3 1 2]);
graze(B.idx.mys,[B.isnek' false false],:,:)=0;

end

%------------------------------
% ZL Diapause Splitting
%------------------------------

if B.diapause
    
% Grazing
    
zlfood = graze(:,B.idx.zl,:,:);
zlloss = graze(B.idx.zl,:,:,:);

zlfood = bsxfun(@times, zlfood, repmat(permute(zlfrac, [4 3 1 2]), [nb+2, 1, 1, 1]));
zlloss = bsxfun(@times, zlloss, repmat(permute(zlfrac, [3 4 1 2]), [1, nb+2, 1, 1]));

graze(:,[B.idx.zl1 B.idx.zl2],:,:) = zlfood;
graze([B.idx.zl1 B.idx.zl2],:,:,:) = zlloss;    

graze(:,B.idx.zl,:,:) = 0;
graze(B.idx.zl,:,:,:) = 0;
    
% Predation
if ~nemflag
    
zlfood = pred(:,B.idx.zl,:,:);
zlloss = pred(B.idx.zl,:,:,:);

zlfood = bsxfun(@times, zlfood, repmat(permute(zlfrac, [4 3 1 2]), [nb+2, 1, 1, 1]));
zlloss = bsxfun(@times, zlloss, repmat(permute(zlfrac, [3 4 1 2]), [1, nb+2, 1, 1]));
    
pred(:,[B.idx.zl1 B.idx.zl2],:,:) = zlfood;
pred([B.idx.zl1 B.idx.zl2],:,:,:) = zlloss;    

pred(:,B.idx.zl,:,:) = 0;
pred(B.idx.zl,:,:,:) = 0;

end
    
end

%------------------------------
% Egestion and excretion
%------------------------------
if nemflag
    
[egest, excrete, graze] =  egeexc(graze, B, nb, nz, nx);

else
    
[egest_in, excrete_in, graze, egest_out, excrete_out] =  egeexc(graze, B, nb, nz, nx, pred_in, pred_out);
   
end

%------------------------------
% Mortality
%------------------------------    

mort = nonpredmort(nemflag, bv.orig, basum1.orig, basum2.orig, bfrac.orig, P, B, G, nb, nz, nx);

%------------------------------
% Ammonification
%------------------------------   

amm = ammonif(bv.orig, P, B, nb, nz, nx);

%------------------------------
% Denitrification
%------------------------------

den = denit(bv.orig, P, B, nb, nz, nx);

%------------------------------
% Nitrification
%------------------------------

nit = nitrif(bv.orig, I(:,:,1), P, B, nb, nz, nx); % I is identical along 3 rd dim.

%------------------------------
% fishing mortality
%------------------------------

fish=fishmort(nemflag, bv.orig, G, B, nb, nz, nx, Arch, it);

% ------------------------------------------------------------------------
%  Flux struct assembly (retain original field names/logic)
%  ------------------------------------------------------------------------

if ~nemflag 
    if B.ecosimppflag
    Flx.amm = amm;
    Flx.den = den;
    Flx.nit = nit;
    Flx.ege_in = egest_in;
    Flx.ege_out = egest_out;
    Flx.exc_in = excrete_in;
    Flx.exc_out = excrete_out;
    Flx.npp = npp;
    Flx.gra = graze;
    Flx.mor = mort;
    Flx.pre_in = pred_in;
    Flx.pre_out = pred_out;
    Flx.fish = fish;
    else
    Flx.amm = amm;
    Flx.den = den;
    Flx.nit = nit;
    Flx.ege_in = egest_in;
    Flx.ege_out = egest_out;
    Flx.exx = exc;
    Flx.exc_in = excrete_in;
    Flx.exc_out = excrete_out;
    Flx.gpp = gpp;
    Flx.gra = graze;
    Flx.mor = mort;
    Flx.pre_in = pred_in;
    Flx.pre_out = pred_out;
    Flx.fish = fish;
    Flx.res = resp;
    end
else
    Flx.amm = amm;
    Flx.den = den;
    Flx.nit = nit;
    Flx.ege = egest;
    Flx.exx = exc;
    Flx.exc = excrete;
    Flx.gpp = gpp;
    Flx.gra = graze;
    Flx.mor = mort;
    Flx.pre_out = pred_out;
    Flx.fish = fish;
    Flx.res = resp;
end   

% Reroute fluxes as indicated by user

if isfield(B, 'reroute') && ~isempty(B.reroute)
    
    for ir = 1:size(B.reroute,1)

        type = B.reroute{ir,1};
        idx1 = B.reroute{ir,2};
        idx2 = B.reroute{ir,3};
        idx3 = B.reroute{ir,4};
        frac = B.reroute{ir,5};

        reroute = Flx.(type)(idx1,idx2,:,:) .* frac;
        Flx.(type)(idx1,idx2,:,:) = Flx.(type)(idx1,idx2,:,:) - reroute;
        Flx.(type)(idx1,idx3,:,:) = Flx.(type)(idx1,idx3,:,:) + reroute;

    end
    
end

%------------------------------------
% Total fluxes and assemble dB/dt
%------------------------------------
if ~nemflag
    
if B.ecosimppflag
    fluxtot_in = Flx.amm + Flx.den + Flx.nit + Flx.ege_in + Flx.exc_in + Flx.npp + ...
              Flx.gra + Flx.mor + Flx.pre_in + Flx.fish;
          
    fluxtot_out = Flx.amm + Flx.den + Flx.nit + Flx.ege_out + Flx.exc_out + Flx.npp + ...
              Flx.gra + Flx.mor + Flx.pre_out + Flx.fish;

else
    fluxtot_in = Flx.amm + Flx.den + Flx.nit + Flx.ege_in + Flx.exx + Flx.exc_in + Flx.gpp + ...
              Flx.gra + Flx.mor + Flx.pre_in + Flx.res + Flx.fish;
          
    fluxtot_out = Flx.amm + Flx.den + Flx.nit + Flx.ege_out + Flx.exx + Flx.exc_out + Flx.gpp + ...
              Flx.gra + Flx.mor + Flx.pre_out + Flx.res + Flx.fish;

end

else
    fluxtot = Flx.amm + Flx.den + Flx.nit + Flx.ege + Flx.exx + Flx.exc + Flx.gpp + ...
              Flx.gra + Flx.mor + Flx.pre_out + Flx.res + Flx.fish;
end

if nemflag
        %nz x nx x nbsv+2 x 1 total flux entering each critter in each box
        fluxin  = permute(sum(fluxtot, 1), [3 4 2 1]);
        %nz x nx x nbsv+2 x 1 total flux leaving each critter in each box
        fluxout = permute(sum(fluxtot, 2), [3 4 1 2]);
else
        %nz x nx x nbsv+2 x 1 total flux entering each critter in each box
        fluxin  = permute(sum(fluxtot_in, 1), [3 4 2 1]);
        %nz x nx x nbsv+2 x 1 total flux leaving each critter in each box
        fluxout = permute(sum(fluxtot_out, 2), [3 4 1 2]);

end


% Final rate of change
db = fluxin(:,:,1:nb,:) - fluxout(:,:,1:nb,:);%nz x nx x nbsv x 1 molN.m-3.s-1

if any(isnan(db(:)))
    warning('BIO:nanInDbdt', 'NaN in dB/dt');
end

% ------------------------------------------------------------------------
%  If diapause: transfer ZL1/ZL2 fluxes back into ZL for reporting
%  ------------------------------------------------------------------------

if B.diapause
    fx = fieldnames(Flx);
    for ii = 1:length(fx)
        Flx.(fx{ii})(B.idx.zl,:,:,:) = Flx.(fx{ii})(B.idx.zl1,:,:,:) + Flx.(fx{ii})(B.idx.zl2,:,:,:);
        Flx.(fx{ii})(:,B.idx.zl,:,:) = Flx.(fx{ii})(:,B.idx.zl1,:,:) + Flx.(fx{ii})(:,B.idx.zl2,:,:);
    end
end

% ------------------------------------------------------------------------
%  Diagnostics (primary production limitation terms)
%  ------------------------------------------------------------------------

Diag.lightlim = Lfc(:,:,[B.idx.ps B.idx.pl]);
Diag.no3lim   = no3lim(:,:,[B.idx.ps B.idx.pl]);
Diag.nh4lim   = nh4lim(:,:,[B.idx.ps B.idx.pl]);
Diag.psmax    = psmax(:,:,[B.idx.ps B.idx.pl]);
Diag.silim    = silim(:,:,B.idx.pl);
Diag.I        = I(:,:,1); % nz x nx 
Diag.kappa    = kappa;
Diag.kp       = kappaP;

order = {'psmax', 'lightlim', 'no3lim', 'nh4lim', 'silim', 'I', ...
         'kappa', 'kp'};
Diag = orderfields(Diag, order);
     

flxcheck = structfun(@(x) any(x(:) < 0), Diag);
if any(flxcheck)
   warning('BIO:negflux', 'Negative flux');
end

% =========================================================================
%  Local helper functions
%  =========================================================================

function q = aydinfrnew(bi, bj, b0, q0, x, d, theta)
% AYDINFRNEW Aydin/EwE functional response (vectorized)
%
% bi, bj : nb x nb arrays (prey and predator biomasses in compatible units)
% b0     : nb x 1 (or 1 x nb) base biomass scaling (Ecopath biomasses)
% q0     : nb x nb base consumption rate (Ecopath rate)
% x, d, theta : shaping parameters (nb x nb)
%
% Returns:
% q : nb x nb consumption (mol N m^-2 s^-1)

nb = size(bi,1);
b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0'; 

yi = bi./b0i; 
yj = bj./b0j; 

q = q0 .* ((x.*yj)./(x - 1 + yj)) .* ((d.*yi.^theta)./(d - 1 + yi.^theta));
q(isnan(q))=0;


function q = ivlev(bi, bj, grmax, kgra, lambda, thresh, temp)
% IVLEV Ivlev grazing functional response with temperature modifier
%
% bi, bj : nb x nb prey/pred arrays (mol N m^-3)
% grmax  : nb x nb maximum grazing
% kgra   : 1 x nb (or nb x 1) temperature sensitivity per predator (broadcasted)
% lambda : 1 x nb (or nb x 1) Ivlev steepness per predator (broadcasted)
% thresh : nb x nb or broadcastable threshold
% temp   : scalar temperature for the box
%
% Returns:
%   q : nb x nb grazing (mol N m^-3 s^-1)

 nb = size(bi,1);
 
 kgra = ones(nb,1) * kgra';
 lambda = ones(nb,1) * lambda';
 
 q = grmax.*exp(kgra.*temp).*(1-exp(lambda.*(thresh-bi))).*bj;
 
