function nit = nitrif(bv, I, P, B, nb, nz, nx)
% NITRIF - Compute nitrification fluxes for the WCVI-E2E model.
%
% This function calculates nitrification (NH₄⁺ → NO₂⁻ or NO₃⁻) based on
% temperature and light conditions, following the model described in:
%   Schmidt M, Eggert A. (2012) A regional 3D coupled ecosystem model of the Benguela upwelling system.
%
% Inputs:
%   bv   - Biomass concentration matrix [nz x nx x nb]
%   I    - Irradiance [nz x nx]
%   P    - Struct with environmental parameters:
%            .T : Temperature [nz x nx]
%   B    - Struct with biological parameters:
%            .Nit0  : Base nitrification rate (scalar)
%            .KNit  : Temperature coefficient (scalar)
%            .I0    : Irradiance threshold (scalar)
%            .KI    : Irradiance half-saturation (scalar)
%            .flux  : Cell array defining fluxes, e.g. {'nit', source, sink}
%   nb   - Number of biological groups
%   nz   - Number of vertical layers
%   nx   - Number of horizontal cells
%
% Output:
%   nit  - Nitrification flux array [nb+2 x nb+2 x nz x nx]
%          Units: typically in molN m⁻³ s⁻¹
%
% Copyright (c) 2026 Virginie Bornarel
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.



% ------------------------------------
% Initialize output array
% ------------------------------------
nit = zeros(nb+2,nb+2,nz,nx);

% ------------------------------------
% Identify nitrification flux pairs
% ------------------------------------
tf = strcmp(B.flux(:,1),'nit');
nitidx = cell2mat(B.flux(tf,[2 3])); 

% ------------------------------------
% Loop through each nitrification reaction
% ------------------------------------
for i=1:size(nitidx,1)
    src = nitidx(i,1);
    snk = nitidx(i,2);

    td = tempdep(B.Nit0, B.KNit, P.T);% nz x nx
    ld = lightdep(I, B.I0, B.KI); % nz x nx

    nit(src,snk,:,:) = reshape(td.*ld.*bv(:,:,src),1,1,nz,nx);

end

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);

% ------------------------
% Light limitation function
% ------------------------

function ld = lightdep(I, I0, KI)
ld = (1 - max(0, (I-I0)./(KI+I-I0)));

