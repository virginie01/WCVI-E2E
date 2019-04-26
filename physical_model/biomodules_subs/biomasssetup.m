function [bv, ba, basum, bfrac, zlfrac, nb, nz, nx] =  biomasssetup(bio, G, B)
%BIOMASSSETUP Biomass calculations for nemurokak/wce modules
%
% [bv, ba, basum, bfrac, zlfrac, nb, nz] =  biomasssetup(bio, A)
%
% The diapause and preyvis options require different processes to "see"
% different fractions of the various functional group biomasses.
%
% Input variables:
%
%   bio:    nz x nx x nbsv array of biomass, values received from ode solver
%
%   G/B:      structure various grid and biological parameters
%
% Output variables:
%
%   bv:     volumetric, per-layer biomass concentration
%           orig:       nz x nx x nb, same as input
%           zlcombo:    nz x nx x nb, ZL1/ZL2 moved to ZL (diapause only)
%           prey:       nb x nb x nz x nx. prey by pred by depth by longitude, 
%                       prey profile for each predator/prey link
%           pred:       nb x nb x nz x nx. prey by pred by depth by longitude, 
%                       predator profile for each predator/prey link 
%
%   ba:     depth-integrated, per-layer biomass, same fields as bv
%
%   basum:  water column-integrated biomass, same fields as bv
%
%   zlfrac: fraction of ZL1 in the ZL total
%
%   nb:     number of biological state variables
%
%   nz:     number of depth layers

% Copyright 2014 Kelly Kearney

%---------------------
% Volumetric biomass
% (per layer, mol/m^3)
%---------------------

% Biomass as passed by mixed layer.  In the diapause case, ZL is 0 and ZL1
% and ZL2 hold the non-diapausing and diapausing copepod populations,
% respectively.

bv.orig = bio; 

% The only alteration here is to move the split ZL1/ZL2 populations into ZL
% if diapause is on, and zero out those subgroups.

bv.zlcombo = bv.orig;
if B.diapause
    bzl = bv.orig(:,:,[B.idx.zl1 B.idx.zl2]);
    
    bv.zlcombo(:,:,B.idx.zl) = sum(bzl,3);
    bv.zlcombo(:,:,[B.idx.zl1 B.idx.zl2]) = 0;
    
    zlfrac = bsxfun(@rdivide, bzl, sum(bzl,3));
    
end

% The predator/prey population varies by diet link, since some predators
% see prey differently.  ZL1/ZL2 is seen as a combined group, and predator
% ability to access prey is set per depth layer.

[nz, nx, nb] = size(bv.orig);

[bv.pred, bv.prey] = deal(zeros(nb,nb,nz,nx));
for iz = 1:nz
    for ix = 1:nx
    btmp = reshape(bv.zlcombo(iz, ix, :),[],1,1);
    
    bv.prey(:,:,iz,ix) = btmp * ones(1,nb);
    bv.pred(:,:,iz,ix) = ones(nb,1) * btmp';
    end
end

% TODO: make sure nekton see other nekton?
    

%---------------------
% Integrated biomass
% (per layer, mol/m^2)
%---------------------

[ba.orig,    basum.orig,    bfrac.orig]    = intoverdepth(bv.orig, G.dz, 1);
[ba.zlcombo, basum.zlcombo, bfrac.zlcombo] = intoverdepth(bv.zlcombo, G.dz, 1);
[ba.pred,    basum.pred,    bfrac.pred]    = intoverdepth(bv.pred, G.dz, 2);
[ba.prey,    basum.prey,    bfrac.prey]    = intoverdepth(bv.prey, G.dz, 2);

%---------------------
% Subfunction: 
% integrate over depth
%---------------------   

function [ba, basum, bfrac] = intoverdepth(bv, dz, flag)
                     
if flag==1    
ba    = bsxfun(@times, bv, dz);                       % nz x nx x nbsv, mol N/m^2 per box
basum = reshape(sum(ba, [1 2]),1,[],1);               % 1 x nbsv, mol N/m^2 total domain
bfrac = bsxfun(@rdivide, ba, reshape(basum,1,1,[]));  % nz x nx x nbsv, fraction of biomass in each box

elseif flag==2
    
dz=permute(dz,[3 4 1 2]);    

ba = bsxfun(@times, bv, dz);                          % nb x nb x nz x nx array, mol N/m^2 per box and per link prey x pred
basum = sum(ba, [3 4]);                               % nb x nb array, mol N/m^2 total domain
bfrac = bsxfun(@rdivide, ba, basum);                  % nb x nb x nz x nx, fraction w/o unit
end
