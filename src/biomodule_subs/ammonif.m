function amm = ammonif(bv, P, B, nb, nz, nx)
% AMMONIF Compute ammonification rates from detrital biomass.
%
% This function calculates the ammonification flux (i.e., mineralization
% of nitrogen from detritus to ammonium) based on oxygen availability
% and temperature using a sigmoidal oxygen dependency and Q10-style 
% temperature response.
%
% Ammonification is computed for all source-sink group pairs for which
% the decomposition rate is defined (vdec > 0).
%
% Inputs:
%   bv   - 3D matrix [nz x nx x nb] of biomass concentration for each group
%   P    - Struct containing environmental fields (P.T = temperature field [nz x nx])
%   B    - Struct containing biological parameters (vdec, Kdec, alphao2, O2)
%   nb   - Number of biological groups
%   nz   - Number of vertical layers
%   nx   - Number of horizontal grid cells
%
% Output:
%   amm  - 4D matrix [nb+2 x nb+2 x nz x nx] of ammonification fluxes
%
% Reference:
%   Schmidt, M., & Eggert, A. (2012). A regional 3D coupled ecosystem model
%   of the Benguela upwelling system. Marine Science Reports, No. 87.
%   https://doi.org/10.12754/MSR-2012-0087
%
% This function was developed for the WCVI-E2E ecosystem model and
% incorporates concepts adapted from decomposition/remineralization
% routines in the original WCE/NEMURO-based framework developed by
% Kelly Kearney.
%
% Modifications and extensions by Virginie Bornarel (2017–2026) include:
%   - adaptation to a 2D coastal upwelling domain
%   - explicit oxygen-dependent remineralization
%   - revised ammonification parameterization
%   - expanded environmental forcing support
%   - updated diagnostics and documentation
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% -------------------------
% Initialize output matrix
% -------------------------
amm = zeros(nb+2,nb+2,nz,nx);

% -------------------------
% Main loop: ammonification
% -------------------------
for isrc = 1:nb
    for isnk = 1:nb
        if B.vdec(isrc,isnk) > 0
            td = tempdep(B.vdec(isrc,isnk), B.Kdec(isrc,isnk), P.T); % nz x nx
            fo2=fdecomp(B.alphao2,B.O2); % nz x nx
            amm(isrc,isnk,:,:) = permute(td .* fo2 .* bv(:,:,isrc), [4 3 1 2]);          
        end
    end
end


%------------------------------------------
% Q10-style temperature dependence function
%------------------------------------------
function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);

% -----------------------------------------------------
% Sigmoidal oxygen-limitation function from Schmidt 2012
% -----------------------------------------------------
function fX = fdecomp (alphax, X)
fX = (1- exp(-2.*alphax.*X))./(1+ exp(-2.*alphax.*X));
