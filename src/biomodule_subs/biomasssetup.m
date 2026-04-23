function [bv, ba1, ba2, basum1, basum2, bfrac, nb, nz, nx, varargout] =  biomasssetup(bio, G, B)
% BIOMASSSETUP Biomass calculations for nemurokak/wce modules
%
% [bv, ba, basum, bfrac, zlfrac, nb, nz, nx] =  biomasssetup(bio, G, B)
%
% The diapause and preyvis options require different processes to "see"
% different fractions of the various functional group biomasses.
%
% Input variables:
%
%   bio:    nz x nx x nbsv array of biomass, values received from ode solver
%
%   G/B:    structure containing various grid and biological parameters
%
% Output variables:
%
%   bv:     volumetric, per-box biomass concentration
%           orig:       nz x nx x nb, same as input
%           zlcombo:    nz x nx x nb, sum of ZL1/ZL2 moved to ZL (diapause only)
%           prey:       nb x nb x nz x nx. prey by pred by depth by longitude, 
%                       prey profile for each predator/prey link
%           pred:       nb x nb x nz x nx. prey by pred by depth by longitude, 
%                       predator profile for each predator/prey link 
%
%   ba:     depth-integrated, per-box biomass, same fields as bv
%
%   basum:  spatially-integrated biomass, same fields as bv
%
%   zlfrac: nz x nx x 2 array. Fraction of ZL1 and ZL2 in the ZL total
%
%   nb:     number of biological state variables
%
%   nz:     number of vertical layers
%
%   nz:     number of horizontal layers

% Copyright 2014 Kelly Kearney

%---------------------
% Volumetric biomass
% (per layer, mol/m^3)
%---------------------

% Biomass as passed by physical model.  In the diapause case, ZL is 0 and ZL1
% and ZL2 hold the non-diapausing and diapausing copepod populations,
% respectively.

bv.orig = bio; 

% The only alteration here is to move the split ZL1/ZL2 populations into ZL
% if diapause is on, and zero out those subgroups.

bv.zlcombo = bv.orig;
if B.diapause
    bzl = bv.orig(:,:,[B.idx.zl1 B.idx.zl2]);
    
    if sum(bzl,3) == 0
    warning('biomasssetup: sum of ZL1 and Zl2 equals to 0')
    end
    
    bv.zlcombo(:,:,B.idx.zl) = sum(bzl,3);
    bv.zlcombo(:,:,[B.idx.zl1 B.idx.zl2]) = 0;
    
    zlfrac = bsxfun(@rdivide, bzl, sum(bzl,3));
    varargout{1} = zlfrac;
end

% 1st version (Kelly's):
% The predator/prey population varies by diet link, since some predators
% see prey differently.  ZL1/ZL2 is seen as a combined group, and predator
% ability to access prey is set per depth layer.
% 2nd version (mine): don't take preyvis into account.

[nz, nx, nb] = size(bv.orig);

[bv.pred, bv.prey] = deal(zeros(nb,nb,nz,nx));
for iz = 1:nz
    for ix = 1:nx
    btmp = reshape(bv.zlcombo(iz, ix, :),[],1,1);
    
    bv.prey(:,:,iz,ix) = btmp * ones(1,nb);
    bv.pred(:,:,iz,ix) = ones(nb,1) * btmp';
    end
end

%---------------------
% Integrated biomass
% (per layer, mol/m^2)
%---------------------

[ba1.orig,   ba2.orig,       basum1.orig,    basum2.orig,    bfrac.orig]    = intoverdepth(bv.orig, G, nz, nx, nb, 1);
[ba1.zlcombo,ba2.zlcombo,    basum1.zlcombo, basum2.zlcombo, bfrac.zlcombo] = intoverdepth(bv.zlcombo, G, nz, nx, nb, 1);
[ba1.pred,   ba2.pred,       basum1.pred,    basum2.pred,    bfrac.pred]    = intoverdepth(bv.pred, G, nz, nx, nb, 2);
[ba1.prey,   ba2.prey,       basum1.prey,    basum2.prey,    bfrac.prey]    = intoverdepth(bv.prey, G, nz, nx, nb, 3);

%---------------------
% Subfunction: 
% integrate over depth
%---------------------   

function [ba1, ba2, basum1, basum2, bfrac] = intoverdepth(bv, G, nz, nx, nb, flag)
                     
if flag==1

%dz1 to calculate true basum coherent w EwE
dz1 = repmat(G.dz,[1,1,nb]); % nz x nx x nbsv, dz in m, modified to accomodate zooplankton densities on slope
dz1(1,2,[1 2 3])=100;dz1(1,2,4)=75;dz1(1,2,5)=50;
dz1(2,2,[1 2 3])=990-100;dz1(2,2,4)=990-75;dz1(2,2,5)=50;
dz1(3,2,5)=0;

%dz2 to distribute predation mortality in each box and coherent with
%concentrations displayed in each box- exlude the extra biomass under the
%MLD on the slope- otherwise the predation mortality in box 22 wouldnt be
%coherent with the concentration displayed in this box
dz2 = repmat(G.dz,[1,1,nb]); % nz x nx x nbsv, dz in m, modified to accomodate zooplankton densities on slope
dz2(2,2,[1 2 3])=990-100;dz2(2,2,4)=990-75;dz2(2,2,5)=50;
%dz2(3,2,5) = 0; I remove it vecause I didn't insert it in versions
%sim1-sim27. don't know why. doesn't change much.

ba1 = bv.*dz1; ba2 = bv.*dz2;                                                              % nz x nx x nbsv, mol N/m^2 per box
areapb = G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2)));                                            % nz x nx, m2 for each box
molpb1 = bsxfun(@times, ba1, areapb); molpb2 = bsxfun(@times, ba2, areapb);                % nz x nx x nbsv, mol N in each box
moltot1 = sum(sum(molpb1,1),2); moltot2 = sum(sum(molpb2,1),2);                            % 1 x 1 x nbsv, total molN
basum1 = reshape((moltot1./G.area),1,[],1); basum2 = reshape((moltot2./G.area),1,[],1);    % 1 x nbsv, mol N/m^2 total domain
bfrac = bsxfun(@rdivide, molpb2, moltot2);                                                 % nz x nx x nbsv, fraction of biomass in each box

elseif flag==2
    
dz1 = zeros(nb,nb,nz,nx);dz2 = zeros(nb,nb,nz,nx);% nb x nb x nz x nx
for i =1:nb
    for j=1:nb
        dz1(i,j,:,:)=G.dz;
        dz2(i,j,:,:)=G.dz;
    end
end
dz1(:,[1 2 3],1,2)=100;dz1(:,4,1,2)=75;dz1(:,5,1,2)=50;
dz1(:,[1 2 3],2,2)=990-100;dz1(:,4,2,2)=990-75;dz1(:,5,2,2)=50;
dz1(:,5,3,2)=0;
dz2(:,[1 2 3],2,2)=990-100;dz2(:,4,2,2)=990-75;dz2(:,5,2,2)=50;%dz2(:,5,3,2)=0;

ba1 = bv.*dz1; ba2 = bv.*dz2;                                              % nb x nb x nz x nx array, mol N/m^2 per box and per link prey x pred
areapb = permute(G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2))),[3 4 1 2]);         % 1 x 1 x nz x nx, m2 for each box
molpb1 = bsxfun(@times, ba1, areapb);molpb2 = bsxfun(@times, ba2, areapb); % nb x nb x nz x nx, mol N in each box
moltot1 = sum(sum(molpb1,3),4);moltot2 = sum(sum(molpb2,3),4);             % nb x nb x 1 x 1, total mol N
basum1 = moltot1./G.area;basum2 = moltot2./G.area;                         % nb x nb array, mol N/m^2 total domain
bfrac = bsxfun(@rdivide, molpb2, moltot2);                                 % nb x nb x nz x nx, fraction w/o unit

elseif flag ==3
    
dz1 = zeros(nb,nb,nz,nx);dz2 = zeros(nb,nb,nz,nx);% nb x nb x nz x nx
for i =1:nb
    for j=1:nb
        dz1(i,j,:,:)=G.dz;
        dz2(i,j,:,:)=G.dz;
    end
end
dz1([1 2 3],:,1,2)=100;dz1(4,:,1,2)=75;dz1(5,:,1,2)=50;
dz1([1 2 3],:,2,2)=990-100;dz1(4,:,2,2)=990-75;dz1(5,:,2,2)=50;
dz1(5,:,3,2)=0;
dz2([1 2 3],:,2,2)=990-100;dz2(4,:,2,2)=990-75;dz2(5,:,2,2)=50;%dz2(5,:,3,2)=0;

ba1 = bv.*dz1; ba2 = bv.*dz2;                                              % nb x nb x nz x nx array, mol N/m^2 per box and per link prey x pred
areapb = permute(G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2))),[3 4 1 2]);         % 1 x 1 x nz x nx, m2 for each box
molpb1 = bsxfun(@times, ba1, areapb);molpb2 = bsxfun(@times, ba2, areapb); % nb x nb x nz x nx, mol N in each box
moltot1 = sum(sum(molpb1,3),4);moltot2 = sum(sum(molpb2,3),4);             % nb x nb x 1 x 1, total mol N
basum1 = moltot1./G.area;basum2 = moltot2./G.area;                         % nb x nb array, mol N/m^2 total domain
bfrac = bsxfun(@rdivide, molpb2, moltot2);                                 % nb x nb x nz x nx, fraction w/o unit


end
  