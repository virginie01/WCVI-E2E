function [gpp, exc, resp, dz, psmax, Lfc, no3lim, nh4lim, silim, I, kappa, kappaP] = primprod(bv, G, P, B, nz, nx, nb)
%PRIMPROD Primary production submodule for wcvie2eode
%
% Computes light attenuation, nutrient limitation (NO3, NH4, SiOH4),
% light limitation, and temperature-dependent gross primary production,
% along with extracellular excretion and respiration fluxes.
%
% INPUTS:
% bv - nz x nx x nb array: current concentrations of state variables
% G - struct: geometric/grid info (G.dz, G.z)
% P - struct: physical/forcing fields (e.g., temperature, PAR)
% B - struct: biological parameters (e.g., Vmax, K values, indices)
% nz, nx - scalars: number of vertical layers and horizontal locations
% nb - scalar: number of state variables
%
% OUTPUTS:
% gpp:    -nb x nb x nz x nx array: gross primary production fluxes (source, sink, z domain, x domain)
%          at current time step
% exc:    -nb x nb x nz x nx array: extracellular excretion fluxes (source, sink, z domain, x domain)
%          at current time step
% resp:   -nb x nb x nz x nx array: phytoplankton respiration fluxes (source, sink, z domain, x domain)
%          at current time step
% dz:     -nz x nx array: depth intervals at current time step
% psmax:  -nz x nx x nb array: max production rate at current temperature
% Lfc:    -nz x nx x nb array: light limitation factor at current time step
% no3lim: -nz x nx x nb array: N03 limitation at current time step
% nh4lim: -nz x nx x nb array: NH4 limitation at current time step
% silim:  -nz x nx x nb array: SiOH4 limitation at current time step
% I:      -nz x nx x nb array: light intensity at current time step
% kappa:  -nz x nx array: total light attenuation coefficient
% kappaP: -nz x nx array: light attenuation due to phytoplankton
%
% This file was derived from the original primprod routine developed
% by Kelly Kearney for the WCE/NEMURO framework and substantially
% modified for the WCVI-E2E coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Major modifications by Virginie Bornarel (2017–2026) include:
%   - adaptation from 1D to 2D light and nutrient limitation calculations
%   - revised grid and forcing interfaces using WCVI-E2E structures
%   - updated light attenuation and phytoplankton self-shading calculations
%   - revised gross primary production, excretion, and respiration bookkeeping
%   - removal of iron-quota terms not used in the WCVI-E2E configuration
%   - expanded documentation and output descriptions
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

%----------------------
% 1. Initial Setup
%----------------------
dz = G.dz;

% Shading phytoplankton biomass (for self-shading effect)
shadingbio = sum(bv(:,:,[B.idx.ps B.idx.pl]), 3); %nz x nx array total phytoplankton concentration molN/m3

%----------------------
% 2. Compute kappaP (biological shading)
%----------------------

nearsurf = [0.5 0.5];    % m, Extend near to surface (since right at surface divides by 0)
zedge = [nearsurf; cumsum(dz)]; %(nz+1) x nx array zedge
pintedge = [shadingbio(1,:).*nearsurf; cumsum(shadingbio).*dz]; %(nz+1) x nx array molN.m-2
kppedge = B.alpha2 .* pintedge./zedge; %(nz+1) x nx array m-1

kappaP=zeros(nz,nx);
for i=1:nx
    kappaP(:,i) = interp1(zedge(:,i), kppedge(:,i), G.z(:,i)); %nz x nx array
end

%----------------------
% 3. Total light attenuation (kappa)
%----------------------
kappa = B.alpha1 + kappaP; %nz x nx array m-1

%----------------------
% 4. Light Intensity
%----------------------
I = zeros(nz,nx);
for i = 1:nx
    I(:,i) = P.par24(i) .* exp(-kappa(:,i) .* G.z(:,i)); % nz x nx array W.m-2
end

%----------------------
% 5. Light Limitation Factor (Lfc)
%----------------------
if B.usesteele
    I = repmat(I,1,1,length(B.Iopt));               %nz * nx * nbsv I is location-dependent 
    Iopt = repmat(ones(nz,nx),1,1,length(B.Iopt));
    Iopt = Iopt.*reshape(B.Iopt,1,1,[]);            %nz x nx x nbsv Iopt is group-dependent
    Lfc = I./Iopt .* exp(1 - I./Iopt);              %nz x nx x nbsv
else
    I = repmat(I,1,1,length(B.alpha));              %nz * nx * nbsv I is location-dependent 
    alpha = repmat(ones(nz,nx),1,1,length(B.alpha));
    alpha = alpha.*reshape(B.alpha,1,1,[]);         %nz x nx x nbsv alpha is group-dependent 
    Vmax = repmat(ones(nz,nx),1,1,length(B.Vmax));
    Vmax = Vmax.*reshape(B.Vmax,1,1,[]);            %nz x nx x nbsv Vmax is group dependent
    Lfc = 1 - exp(-alpha .* I./Vmax);               %nz x nx x nbsv
end

%----------------------
% 6. Initialize Outputs
%----------------------
[no3lim, nh4lim, silim, fratio, psmax] = deal(zeros(nz, nx, nb));
[gpp, exc, resp] = deal(zeros(nb+2, nb+2, nz, nx));

%----------------------
% 7. Loop over Phytoplankton Groups
%----------------------
for ib = [B.idx.ps B.idx.pl]
    
    % Macronutrient limitation
    
    no3lim(:,:,ib) = mich(bv(:,:,B.idx.no3), B.Kno3(ib)) .* exp(-B.pusai(ib) .* bv(:,:,B.idx.nh4));
    nh4lim(:,:,ib) = mich(bv(:,:,B.idx.nh4), B.Knh4(ib));
    silim(:,:,ib)  = mich(bv(:,:,B.idx.sioh4), B.Ksi(ib)); % Always 1 if Ksi == 0, i.e. not si-dependant
    
    % F-ratio (fraction new production)   
    
    fratio(:,:,ib) = no3lim(:,:,ib)./(no3lim(:,:,ib) + nh4lim(:,:,ib));% nz x nx array
      
    % Overal nutrient limitation
    
    % nutlim = min([(no3lim(:,ib) + nh4lim(:,ib)) silim(:,ib) felim(:,ib)], [], 2);
    if ib == B.idx.ps
        nutlim = no3lim(:,:,ib) + nh4lim(:,:,ib); %nz x nx array
    else   
        nutlim = min(cat(3,(no3lim(:,:,ib) + nh4lim(:,:,ib)), silim(:,:,ib)), [], 3); %nz x nx array
    end
    
    % Maximum possible production rate at current temps
    
    psmax(:,:,ib) = tempdep(B.Vmax(ib), B.Kgpp(ib), P.T); %nz x nx array
    
    % Nitrogen uptake
    
    totnuptake = psmax(:,:,ib) .* Lfc(:,:,ib) .* nutlim .* bv(:,:,ib);%nz x nx array
    
    no3uptake = totnuptake .* fratio(:,:,ib); %nz x nx array
    nh4uptake = totnuptake .* (1 - fratio(:,:,ib)); %nz x nx array
    
    % Uptake of SiOH4
    
    if ib == B.idx.pl
        siuptake = (no3uptake + nh4uptake) .* B.RSiN; %nz x nx array
    else
        siuptake = zeros(nz,nx); %
    end
    
    % Gross primary production fluxes
    
    gpp(B.idx.no3,ib,:,:) = no3uptake;
    gpp(B.idx.nh4,ib,:,:) = nh4uptake;
    if ib == B.idx.pl
        gpp(B.idx.sioh4,B.idx.plsi,:,:) = siuptake;
    end

    % Extracellular excretion (N) is proportional to uptake
    
    excn = B.gamma(ib) .* (no3uptake + nh4uptake);
    exc(ib,B.idx.don,:,:) = excn;
        
    % Phytoplankton respiration (N)
    
    respn = tempdep(B.res0(ib), B.Kres(ib), P.T) .* bv(:,:,ib); %nz x nx array
    
    resp(ib,B.idx.no3,:,:) = respn .* fratio(:,:,ib);
    resp(ib,B.idx.nh4,:,:) = respn .* (1 - fratio(:,:,ib));
    
    % Excretion and respiration of silica proportional to N
    
    if ib == B.idx.pl
        excsi = excn * B.RSiN;
        exc(B.idx.plsi,B.idx.sioh4,:,:) = excsi;
        
        respsi = respn * B.RSiN;
        resp(B.idx.plsi, B.idx.sioh4,:,:) = respsi;
    end
    
    
end

%------------------------
% Michaelis-Menten uptake
% limitation
%------------------------

function lim = mich(x, kx)
lim = x./(x + kx);
lim(isnan(lim)) = 0; % if x == 0 and kx == 0

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);
