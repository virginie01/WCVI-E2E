function [db, Flx, Diag] = wceode(time, bio, G, P, B)
%WCEODE Water column ecosystem model, main ODE function
%
% [db, Flx, Diag] = wceode(time, bio, A)
%
% Source/sink ODE function for wce module.  See biomodules/wce for details;
% this function is designed to be called by ODE solvers.


bio = max(bio, 0); % Negative biomass treated as 0

%------------------------------
% Set up various bio arrays
%------------------------------


% [bv,    bvin,    bvpre,    bvgra, ...
%  basum, basumin, basumpre, basumgra, ...
%  bfrac, bfracin, bfracpre, bfracgra, ...
%         zlfracin, zlfracpre, zlfracgra, nz, nb] = biomasssetup(bio, A);
    
[bv, ba, basum, bfrac, zlfrac, nb, nz, nx] =  biomasssetup(bio, G, B);

%------------------------------
% Photosynthesis-related
% fluxes (gpp, exc, resp)
%------------------------------

[gpp, exc, resp, dz, psmax, Lfc, no3lim, nh4lim, silim, I, ...
    kappa, kappaP] = primprod(bv.orig, G, P, B, nz, nx, nb);


%------------------------------
% Aydin Ecosim primary 
% production
%------------------------------

if B.ecosimppflag  
    npp = zeros(nb+2,nb+2,nz,nx);
    istoodeep = G.z < -100;
    for ib = [B.idx.ps B.idx.pl]
        h = exp(-2) + 1;  % TODO
        ybio = bv.orig(:,:,ib) ./ B.b0v(ib);
        npp(B.idx.mys,ib,:,:) = B.p0v(ib) .* (h.*ybio)./(h-1+ybio);
        npp(B.idx.mys,ib,istoodeep,:) = 0;
    end
end
   
%------------------------------
% Predation (involves nekton)
%------------------------------

pred = zeros(nb+2,nb+2,nz,nx);

% Predation flux based on aydin ecosim functional response
% (ZL never predate, so don't have to worry about separating in/out here)

for iz=1:nz
    for ix=1:nx
    pred(1:nb, 1:nb, iz, ix) = aydinfrnew(bv.prey(:,:,iz,ix), bv.pred(:,:,iz,ix),...
    B.b0v, B.q0v, B.x, B.d, B.theta); % nb x nb, mol N m^-3 s^-1
    end
end

% Involves only nekton

for iprey = 1:nb
    for ipred = 1:nb
        if B.links(iprey,ipred) ~= 2 && B.links(iprey,ipred) ~= 3 &&...
               B.links(iprey,ipred) ~= 5 && B.links(iprey,ipred) ~= 7 
           pred(iprey,ipred,:) = 0;
        end
    end
end

%------------------------------
% Grazing (involves only 
% plankton, and includes temp
% influence)
%------------------------------

% grazing flux based on aydin ecosim functional response with temperature
% dependence

graze = zeros(nb+2,nb+2,nz,nx);
for iz = 1:nz
    for ix= 1:nx
    graze(1:nb,1:nb,iz,ix) = aydinfrtempnew(bv.prey(:,:,iz,ix), bv.pred(:,:,iz,ix),...
        B.b0v, B.q0vat0, B.x, B.d, B.theta, B.Kgra, P.T(iz,ix)); % mol N m^-3 s^-1
    end    
end

% Involves only plankton

for iprey = 1:nb
    for ipred = 1:nb
        if B.links(iprey,ipred) ~= 1 && B.links(iprey,ipred) ~= 4 &&...
               B.links(iprey,ipred) ~= 6 
           graze(iprey,ipred,:,:) = 0;
        end
    end
end

% If diapause flag is on, split the ZL intake and loss proportionally
% across the ZL1 and ZL2 groups

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
    
zlfood = pred(:,B.idx.zl,:,:);
zlloss = pred(B.idx.zl,:,:,:);

zlfood = bsxfun(@times, zlfood, repmat(permute(zlfrac, [4 3 1 2]), [nb+2, 1, 1, 1]));
zlloss = bsxfun(@times, zlloss, repmat(permute(zlfrac, [3 4 1 2]), [1, nb+2, 1, 1]));
    
pred(:,[B.idx.zl1 B.idx.zl2],:,:) = zlfood;
pred([B.idx.zl1 B.idx.zl2],:,:,:) = zlloss;    

pred(:,B.idx.zl,:,:) = 0;
pred(B.idx.zl,:,:,:) = 0;
    
end

%------------------------------
% Egestion and excretion
%------------------------------
    
[egest, excrete, graze] =  egeexc(pred, graze, B, nb, nz, nx);

%------------------------------
% Mortality
%------------------------------    

mort = nonpredmort(false, bv.orig, P, B, nb, nz, nx);

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

nit = nitrif(bv.orig, G, P, B, nb, nz, nx);

%------------------------------
% fishing mortality
%------------------------------

fmort=fishmort(false, bv.orig, bfrac.orig, G, B, nb, nz, nx);

%------------------------------
% Reroute
%------------------------------ 

% Gather together flux terms

if B.ecosimppflag
    Flx.amm = amm;
    Flx.den = den;
    Flx.nit = nit;
    Flx.ege = egest;
    Flx.exc = excrete;
    Flx.npp = npp;
    Flx.gra = graze;
    Flx.mor = mort;
    Flx.pre = pred;
    Flx.fmort = fmort;

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
    Flx.pre = pred;
    Flx.fmort = fmort;
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

%------------------------------
% Total fluxes
%------------------------------  

if B.ecosimppflag
    fluxtot = Flx.amm + + Flx.den + + Flx.nit + Flx.ege + Flx.exc + Flx.npp + ...
              Flx.gra + Flx.mor + Flx.pre + Flx.fmort;
else
    fluxtot = Flx.amm + Flx.den + Flx.nit + Flx.ege + Flx.exx + Flx.exc + Flx.gpp + ...
              Flx.gra + Flx.mor + Flx.pre + Flx.res + Flx.fmort;
end

%nz x nx x nbsv total flux entering each critter in each box
fluxin  = permute(sum(fluxtot, 1), [3 4 2 1]);
%nz x nx x nbsv total flux leaving each critter in each box
fluxout = permute(sum(fluxtot, 2), [3 4 1 2]);

% Final rate of change
db = fluxin(:,:,1:nb) - fluxout(:,:,1:nb);%nz x nx x nbsv

if any(isnan(db(:)))
    warning('WCE:nanInDbdt', 'NaN in dB/dt');
end

% If diapause flag is set, need to transfer the fluxes from the ZL1/ZL2
% groups to ZL

if B.diapause
    fx = fieldnames(Flx);
    for ii = 1:length(fx)
        Flx.(fx{ii})(B.idx.zl,:,:,:) = Flx.(fx{ii})(B.idx.zl1,:,:,:) + Flx.(fx{ii})(B.idx.zl2,:,:,:);
        Flx.(fx{ii})(:,B.idx.zl,:,:) = Flx.(fx{ii})(:,B.idx.zl1,:,:) + Flx.(fx{ii})(:,B.idx.zl2,:,:);
    end
end

% Diagnostics (all intermediate fluxes)

Diag.lightlim = Lfc(:,:,[B.idx.ps B.idx.pl]);
Diag.no3lim   = no3lim(:,:,[B.idx.ps B.idx.pl]);
Diag.nh4lim   = nh4lim(:,:,[B.idx.ps B.idx.pl]);
Diag.psmax    = psmax(:,:,[B.idx.ps B.idx.pl]);
Diag.silim    = silim(:,:,B.idx.pl);
Diag.I        = I;
Diag.kappa    = kappa;
Diag.kp       = kappaP;

order = {'psmax', 'lightlim', 'no3lim', 'nh4lim', 'silim', 'I', ...
         'kappa', 'kp'};
Diag = orderfields(Diag, order);
     

flxcheck = structfun(@(x) any(x(:) < 0), Diag);
if any(flxcheck)
   warning('WCE:negflux', 'Negative flux');
end

%------------------------------
% Aydin functional response
%------------------------------

function q = aydinfrnew(bi, bj, b0, q0, x, d, theta)

nb = size(bi,1);
b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0'; 

yi = bi./b0i; 
yj = bj./b0j; 

q = q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);


function q = aydinfrtempnew(bi, bj, b0, q0, x, d, theta, kgra, temp)

nb = size(bi,1);

b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0';

kgra = ones(nb,1) * kgra';

yi = bi./b0i;
yj = bj./b0j;

q = exp(kgra.*temp) .* q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);
