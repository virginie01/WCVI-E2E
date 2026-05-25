function den = denit(bv, P, B, nb, nz, nx)
% DENIT Compute denitrification flux for the WCE model.
%
% This function implements the denitrification process based on the 
% coupled ecosystem model described in:
% Schmidt M, Eggert A. (2012). A regional 3D coupled ecosystem model of the
% Benguela upwelling system. Marine Science Reports No. 87.
%
% Input:
%   bv   - Biomass concentrations [nz x nx x nb]
%   P    - Struct with environmental parameters:
%            .T : Temperature [nz x nx]
%   B    - Struct with biological parameters:
%            .vdec     : Base decay rate matrix [nb x nb]
%            .Kdec     : Temperature sensitivity matrix [nb x nb]
%            .alphao2  : Oxygen half-saturation parameter (scalar)
%            .alphaNo3 : Nitrate half-saturation parameter (scalar)
%            .alphaNH4 : Ammonium inhibition parameter (scalar)
%            .O2       : Oxygen concentration [nz x nx]
%            .idx.no3  : Index of NO3 group in bv
%            .idx.nh4  : Index of NH4 group in bv
%            .flux     : Cell array of flux rules: {'name', src, sink}
%   nb   - Number of biological groups
%   nz   - Number of vertical layers
%   nx   - Number of horizontal boxes
%
% Output:
%   den  - Denitrification flux array [nb+2 x nb+2 x nz x nx]
%
% This function was developed for the WCVI-E2E ecosystem model and
% incorporates concepts adapted from decomposition/remineralization
% routines in the original WCE/NEMURO framework developed by
% Kelly Kearney.
%
% Modifications and extensions by Virginie Bornarel (2017–2026) include:
%   - implementation of oxygen-dependent denitrification
%   - nitrate and ammonium limitation terms
%   - adaptation to a 2D coastal upwelling domain
%   - revised environmental forcing support
%   - expanded process-specific flux bookkeeping
%   - updated diagnostics and documentation
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% ------------------------------------
% Initialize output array
% ------------------------------------

den = zeros(nb+2,nb+2,nz,nx);

% ------------------------------------
% Environmental limitation factors
% ------------------------------------

fo2 = fdecomp(B.alphao2,B.O2);                 % [nz x nx] - Oxygen inhibition
fNo3 = fdecomp(B.alphaNo3, bv(:,:,B.idx.no3)); % [nz x nx] - Nitrate availability
fNH3 = fdecomp(B.alphaNH4, bv(:,:,B.idx.nh4)); % [nz x nx] - Ammonium inhibition

% ------------------------------------
% Locate denitrification-related fluxes
% ------------------------------------
tf = strcmp(B.flux(:,1),'den');      % Find "den" fluxes
denidx = cell2mat(B.flux(tf,[2 3])); % Extract [src sink] pairs

% ------------------------------------
% Loop through denitrification fluxes
% ------------------------------------
for i=1:size(denidx,1)
    src = denidx(i,1); % Donor group
    snk = denidx(i,2); % Recipient group

    % Temperature effect (Q10-style)
    td = tempdep(B.vdec(src,snk), B.Kdec(src,snk), P.T); % [nz x nx]
    % Store flux in output array
    den(src,snk,:,:) = permute(td.*(1-fo2).*fNo3.*(1-fNH3).*bv(:,:,src), [4 3 1 2]); % [src x snk x nz x nx]

end


%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);

% -------------------------------------
% Oxygen/nutrient availability limitation
% Based on: Schmidt & Eggert (2012)
% -------------------------------------

function fX = fdecomp (alphax, X)
fX = (1- exp(-2.*alphax.*X))./(1+ exp(-2.*alphax.*X));

