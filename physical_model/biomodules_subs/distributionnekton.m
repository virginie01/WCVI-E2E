function nekpostmov = distributionnekton(nekpremov, nekdist, G, nz, nx)

% NEKPOSTMOV enables to distribute nekton groups according to the input 
% "nekdist" at the beginning of each time step.
%
% Input variables:
%
%     nekpremov: nz x nx x nnek array. Concentration of each nekton group
%                in mol N.m-3 at next time step before any movement/redistribution
%                occurs.
%
%     nekdist:   nnek x 12 cell array. Each cell is a nz x nx matrix holding 
%                percentage of total biomass in each spatial box for a given 
%                nekton group (row) and month(column). Nekton groups are
%                organized in the same order as in Ecopath.
%
%
%     G:         Structure holding grid parameters
%
%
% Output variable
% 
%     nekpostmov: nz x nx x nnek array. Concentration of each nekton group
%                 in mol N.m-3 at next time step after spatial redistribution 
%                 occurs. 
%

nnek = size(nekdist,1);
nekpostmov = zeros(nz,nx,nnek);


% transform molN.m-3 into molN
mol_nek = nekpremov .* (G.dz .* G.dx .* (G.area./sum(G.dx(1,:)))); %nz x nx x nnek
% sum moles across the domain for each nekton group
sum_mol_nek = permute(sum(mol_nek, [1 2]), [1 3 2]); % 1 x nnek

% find month corresponding to next time step

daten = (datenum(G.sdate) + G.t./86400) + G.dt./86400;

datev = datevec(daten);
month = datev(2);

% redistribute nekton biomasses according to time of the year and "nekdist"
% input

for inek = 1:size(nekdist,1)
    nekpostmov(:,:,inek) = sum_mol_nek(inek).*nekdist{inek,month}; %nz x nx molN
    nekpostmov(:,:,inek) = nekpostmov(:,:,inek)./(G.dz .* G.dx .* (G.area./sum(G.dx(1,:))));%nz x nx molN.m-3
end

