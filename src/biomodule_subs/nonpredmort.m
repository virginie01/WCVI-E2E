function mort = nonpredmort(nemflag, bv, basum1, basum2, bfrac, P, B, G, nb, nz, nx)
% NONPREDMORT Compute natural (non-predatory) mortality for plankton and nekton.
%
% SYNTAX:
%   mort = nonpredmort(nemflag, bv, basum1, basum2, bfrac, P, B, G, nb, nz, nx)
%
% INPUTS:
%   nemflag   - Logical flag for NEMURO-based formulation
%   bv        - Biomass matrix [nz x nx x nb]
%   basum1    - Group-specific vertically integrated biomass (used for nekton mortality) [1 x nb]
%   basum2    - Group-specific vertically integrated biomass (used for plankton mortality) [1 x nb]
%   bfrac     - 3D matrix with biomass fraction per box [nz x nx x nb]
%   P         - Struct with environmental variables, e.g., temperature P.T [nz x nx]
%   B         - Biological parameters struct (contains m0coef, m0exp, mor0, Kmor, isphy, iszoo, isnek, idx)
%   G         - Grid structure (area, dz, dx)
%   nb        - Number of biological groups
%   nz, nx    - Number of vertical and horizontal layers
%
% OUTPUT:
%   mort      - Mortality flux matrix [nb+2 x nb+2 x nz x nx], from group to PON/Opal
%
% REMARKS:
%   - If `nemflag` is true, quadratic temperature-dependent mortality is used.
%   - Otherwise, a flexible quadratic mortality function is used with biomass integration and volume scaling.
%
% This file was initially derived from the original nonpredmort
% routine developed by Kelly Kearney for the WCE/NEMURO framework
% and subsequently extended for the WCVI-E2E coastal upwelling
% ecosystem model.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Major modifications and extensions by Virginie Bornarel (2017–2026)
% include:
%   - adaptation from 1D to 2D spatial mortality calculations
%   - separate mortality treatment for plankton and nekton groups
%   - biomass integration and redistribution across WCVI grid cells
%   - revised volume-scaling and domain geometry handling
%   - expanded mortality bookkeeping and diagnostics
%   - updated documentation and process representation
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.


mort = zeros(nb+2,nb+2,nz,nx);

% N

if nemflag
 
    % Temperature-dependent quadratic mortality (NEMURO style)
    mor0=repmat(reshape(B.mor0,1,1,[]),nz,nx,1); % nz x nx x nbsv group-dependent
    temp=repmat(P.T,1,1,length(B.mor0)); % nz x nx x nbsv location-dependent
    kmor=repmat(reshape(B.Kmor,1,1,[]),nz,nx,1); % nz x nx x nbsv group-dependent 
    mortn = bsxfun(@times, bv.^2, tempdep(mor0, kmor, temp)); % nz x nx x nbsv
    mort(1:nb,B.idx.pon,:,:) = permute(mortn, [3 4 1 2]);

else

    % General quadratic mortality
    mortn = zeros(nz,nx,nb);

    % Nektonic mortality (non-spatialized)
    mortnek = bsxfun(@times,bsxfun(@power,basum1,reshape(B.m0exp,1,[])),reshape(B.m0coef,1,[])); % 1 x nb mortality mol N.m-2.s-1
    mortnek = mortnek.*G.area;% 1 x nb mortality in mol N.s-1
    isnek = reshape(B.isnek,1,1,[]);
    mortn(3,1,isnek) = reshape(mortnek(B.isnek')./(G.dz(3,1).*G.dx(3,1).*(G.area./(G.dx(1,1)+G.dx(1,2)))),1,1,[]); % nz x nx x nb molN.m-3.s-1

    mortplankton = bsxfun(@times,bsxfun(@power,basum2,reshape(B.m0exp,1,[])),reshape(B.m0coef,1,[])); % 1 x nb mortality mol N.m-2.s-1
    mortplankton = mortplankton.*G.area;% 1 x nb mortality in mol N.s-1
    mortplankton = bsxfun(@times, bfrac, reshape(mortplankton,1,1,[])); %nz x nx x nb mortality in molN.s-1 in each box

    % Volume per box for nekton
    vpb = G.dz.*G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2))); %nz x nx m3
    vpb = repmat(vpb,1,1,nb);%nz x nx x nb volume per box m3
    vpb(2,2,[1 2 3]) = (990-100).*G.dx(2,2).*(G.area./(G.dx(1,1)+G.dx(1,2)));
    vpb(2,2,4) = (990-75).*G.dx(2,2).*(G.area./(G.dx(1,1)+G.dx(1,2)));
    vpb(2,2,5) = 50.*G.dx(2,2).*(G.area./(G.dx(1,1)+G.dx(1,2)));
    %vpb(3,2,5) = 0; to remove to avoid NaN thereafter

    % Planktonic mortality (spatialized)
    mortplankton = mortplankton./vpb; % nz x nx x nb mortality in molN.m-3.s-1

    isplankton = B.isphy | B.iszoo;
    isplankton = reshape(isplankton,1,1,[]);

    mortn(:,:,isplankton) = mortplankton(:,:,isplankton);% nz x nx x nb mortality molN.m-3.s-1

    mort(1:nb,B.idx.pon,:,:) = permute(mortn, [3 4 1 2]);

end

% Si

mort(B.idx.plsi,B.idx.opal,:,:) = permute(mortn(:,:,B.idx.pl) .* B.RSiN, [4 3 1 2]);

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);