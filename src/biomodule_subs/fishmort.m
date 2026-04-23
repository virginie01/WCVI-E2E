function fish=fishmort(nemflag, bv, G, B, nb, nz, nx, Arch, it)

% forcing time-series data can be either forcing fishing mortality by pool
% (4)or forced catches (-6)

fish = zeros(nb+2,nb+2,nz,nx);

if ~nemflag

ndata=size(B.fish,2);

for data = 1:ndata
    poolidx = B.fish{1,data};
    type = B.fish{2,data};
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
F = -log(1-(c./biomass));
F(isnan(F)) = 0;  