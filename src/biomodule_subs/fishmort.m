function fish=fishmort(nemflag, bv, G, B, nb, nz, nx, Arch, it)
%FISHMORT Apply fishing mortality or catch forcing to biomass pools
%
% This function handles fishing mortality depending on whether fishing is 
% prescribed as instantaneous mortality rates (type 4) or as total catches 
% (type -6), which are converted to mortality rates.
%
% Inputs:
%   nemflag - Boolean, true if NEMURO standalone model is used
%   bv      - Biomass matrix [nz x nx x nb]
%   G       - Grid struct with .dz, .dx, .area, .dt, .t
%   B       - Biological parameters including .fish and .idx.fish
%   nb      - Number of biological groups
%   nz, nx  - Vertical and horizontal dimensions
%   Arch    - Struct array for archiving diagnostic outputs
%   it      - Current simulation time index
%
% Output:
%   fish    - Fishing mortality fluxes [nb+2 x nb+2 x nz x nx]
%
% Copyright (c) 2026 Virginie Bornarel
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% Initialize output
fish = zeros(nb+2,nb+2,nz,nx);

if ~nemflag

ndata=size(B.fish,2); % Number of fishing forcing entries

for data = 1:ndata
    poolidx = B.fish{1,data};
    type = B.fish{2,data};
    % Case 1: Forced catch (in tonnes, converted to molN)
           if type == -6
              Ctmp = B.fish{3,data};%catch (in t) at current time step
              Ctmp = Ctmp .*1885; % transform catch in mole N
               
              % check whether total catch < total biomass
              Biomass = bv(3,1,poolidx).*(G.dz(3,1).*G.dx(3,1).*(G.area./sum(G.dx(1,:))));
              if gt(Ctmp,Biomass)
                  for io = 1:length(Arch)
                      isdl = cellfun(@(x) isequal([nz nx], size(x)), Arch(io).avg);
                      iss = cellfun(@(x) isscalar(x), Arch(io).avg);
    
                      start = cell(size(Arch(io).avg));
                      [start{isdl}] = deal([1 1 Arch(io).bin(it)]);
                      [start{iss}] = deal(Arch(io).bin(it));
    

                      for iv = 1:length(Arch(io).avg)
                        ncwrite(Arch(io).file, Arch(io).datatoarchive{iv,2}, Arch(io).avg{iv}, start{iv});
                      end
                  end
                  
                  error('Total catches are higher than total available biomass. \nGroup: %d, time: %d s since simulation start time',...
                  poolidx, G.t)
                  
              end
               
               % distribute overall catch in box31 given that all nekton
               % groups are accumulated in this box
               C = Ctmp; % (molN)
               C = C./(G.dz(3,1).*G.dx(3,1).*(G.area./sum(G.dx(1,:)))); %molN.m-3 
               C = C./G.dt; %molN.m-3.s-1 
               
               F = fishingmortality(C, bv(3,1,poolidx)); % s-1
               fish(poolidx,B.idx.fish,3,1) = F.* bv(3,1,poolidx); % molN.m-3.s-1
               
           elseif type == 4 %fishing mortality
               F=B.fish{3,data};%s-1 scalar
               fish(poolidx,B.idx.fish,3,1) = F.*bv(3,1,poolidx);
    
           end
    
end
end

function F = fishingmortality(c, biomass)
%FISHINGMORTALITY Compute instantaneous fishing mortality rate from catch
%
% Inputs:
%   c       - Catch rate [molN m⁻³ s⁻¹]
%   biomass - Biomass [molN m⁻³]
%
% Output:
%   F       - Instantaneous mortality rate [s⁻¹]
F = -log(1-(c./biomass));
F(isnan(F)) = 0;  