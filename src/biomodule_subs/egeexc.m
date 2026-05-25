function [egest_in, excrete_in, graze, varargout] = egeexc(graze, B, nb, nz, nx, varargin)
%EGEXC Computes egestion and excretion fluxes for zooplankton and nekton groups.
%
% SYNTAX:
%   [egest_in, excrete_in, graze] = egeexc(graze, B, nb, nz, nx)
%   [egest_in, excrete_in, graze, egest_out, excrete_out] = egeexc(graze, B, nb, nz, nx, pred_in, pred_out)
%
% DESCRIPTION:
%   Computes spatialized egestion and excretion fluxes from predator-prey
%   interactions. Handles both plankton-only grazing and plankton+nekton
%   predation if `pred_in` and `pred_out` are provided.
%
% INPUTS:
%   graze      - [nb+2 x nb+2 x nz x nx] Grazing fluxes (predator eats prey)
%   B          - Structure with biological parameters (e.g., gs, ge, alphaeg, beta, idx, RSiN)
%   nb         - Number of biological groups
%   nz         - Number of depth layers
%   nx         - Number of horizontal grid cells
%
% OPTIONAL INPUTS:
%   pred_in    - [nb+2 x nb+2 x nz x nx] Predation fluxes (non-spatialized), sink side
%   pred_out   - [nb+2 x nb+2 x nz x nx] Predation fluxes (spatialized when prey is P or Z), source side
%
% OUTPUTS:
%   egest_in       - [nb+2 x nb+2 x nz x nx] Egested material received by detritus (PON/Opal)
%   excrete_in     - [nb+2 x nb+2 x nz x nx] Excreted material received by NH4
%   graze          - [nb+2 x nb+2 x nz x nx] Original grazing flux (same as input)
%
% OPTIONAL OUTPUTS (only if predation fluxes provided):
%   egest_out      - [nb+2 x nb+2 x nz x nx] Egested material leaving predators
%   excrete_out    - [nb+2 x nb+2 x nz x nx] Excreted material leaving predators
%
% This file was derived from the original egeexc routine developed
% by Kelly Kearney for the WCE/NEMURO framework and subsequently
% extended for the WCVI-E2E coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Modifications and extensions by Virginie Bornarel (2017–2026) include:
%   - adaptation from a 1D to 2D spatial framework
%   - explicit handling of spatialized predation fluxes
%   - separate inflow/outflow flux bookkeeping
%   - revised nekton-plankton coupling
%   - expanded silica and detrital flux handling
%   - modularized optional predation support
%   - updated diagnostics and documentation
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.


% ---------------------------
% 1. Initialize flux totals
% ---------------------------
if (length(varargin) == 2)
    pred_in = varargin{1}; %nb+2 x nb+2 x nz x nx- nekton side (i.e. sink), non spatialized, everything jammed in box 31
    pred_out = varargin{2}; %nb+2 x nb+2 x nz x nx - nekton+plankton side (i.e. source), partly spatialized, plankton losses spatially distributed and nekton losses jammed in box 31
    alleat_in = pred_in + graze; % (never overlap) nb+2 x nb+2 x nz x nx, N-Z fluxes non spatialized
    alleat_out = pred_out + graze; % (never overlap) nb+2 x nb+2 x nz x nx, N-Z fluxes spatialized
    
else
    alleat = graze;
end

if (length(varargin) == 2)
    
    alleat_in = permute(sum(alleat_in,1), [3 4 2 1]);  %  nz x nx x nb+2, total eaten by each predator- all fluxes involving nekton groups are contained in box 31
    alleat_out = permute(sum(alleat_out,1), [3 4 2 1]); %  nz x nx x nb+2, total eaten by each predator- fluxes involving nekton groups are partially spatialized, N-Z are spatialized and N-N in box 31

else
    
    alleat = permute(sum(alleat,1), [3 4 2 1]); %  nz x nx x nb+2, total eaten by each predator

end

[egest_in, egest_out, excrete_in, excrete_out] = deal(zeros(nb+2,nb+2,nz,nx));

% ---------------------------
% 2. Egestion & Excretion
% ---------------------------

if (length(varargin) == 2)
    egetmp_in = bsxfun(@times, alleat_out, reshape([B.gs' 0 0],1,1,[])); % nz x nx x nb+2, egestion that enters PON- Z-Z and N-Z fluxes spatialized
    egetmp_out = bsxfun(@times, alleat_in, reshape([B.gs' 0 0],1,1,[])); % nz x nx x nb+2, egestion that leaves living critters- spatialized for Z and in box 31 for N
    exctmp_in = bsxfun(@times, alleat_out, reshape([(1 - B.gs' - B.ge') 0 0],1,1,[])); % nz x nx x nb+2 excretion that enters NH4 - Z-Z and N-Z fluxes spatialized
    exctmp_out = bsxfun(@times, alleat_in, reshape([(1 - B.gs' - B.ge') 0 0],1,1,[])); % nz x nx x nb+2 excretion that leaves living critters - spatialized for Z and in box 31 for N
else
    egetmp = bsxfun(@times, alleat, reshape([(1-B.alphaeg') 0 0],1,1,[])); % nz x nx x nb+2, egestion
    exctmp = bsxfun(@times, alleat, reshape([(B.alphaeg' - B.beta') 0 0],1,1,[])); % nz x nx x nb+2 excretion
end

if (length(varargin) == 2)
    egest_in(:,B.idx.pon,:,:) = permute(egetmp_in,[3 4 1 2]);
    egest_out(:,B.idx.pon,:,:) = permute(egetmp_out,[3 4 1 2]);
    excrete_in(:,B.idx.nh4,:,:) = permute(exctmp_in,[3 4 1 2]);
    excrete_out(:,B.idx.nh4,:,:) = permute(exctmp_out,[3 4 1 2]);
else
    egest_in(:,B.idx.pon,:,:) = permute(egetmp,[3 4 1 2]);
    excrete_in(:,B.idx.nh4,:,:) = permute(exctmp,[3 4 1 2]);
end
  
% ---------------------------
% 3. Silicate Egestion (opal)
% ---------------------------
  
% Assume no critters can assimilate Si, all egested immediately to opal

if (length(varargin) == 2)
    grasi = graze(B.idx.pl,:,:,:) + pred_out(B.idx.pl,:,:,:); % 1 x nb+2 x nz x nx
else
    grasi = graze(B.idx.pl,:,:,:);
end

grasi = permute(grasi, [3 4 2 1]);                    % nz x nx x nb+2 
grasi = sum(grasi,3) .* B.RSiN;                       % nz x nx total amount of plSi eaten in each box

egest_in(B.idx.plsi,B.idx.opal,:,:) = reshape(grasi,1,1,nz,nx);
egest_out(B.idx.plsi,B.idx.opal,:,:) = reshape(grasi,1,1,nz,nx);

% ---------------------------
% 4. Optional Outputs
% ---------------------------

if (length(varargin) == 2)
    varargout{1} = egest_out;
    varargout{2} = excrete_out;
end

