function fmort=fishmort(nemflag, bv, bfrac, G, B, nb, nz, nx)

% forcing time-series data can be either forcing fishing mortality by pool
% (4)or forced catches (-6)

fmort = zeros(nb+2,nb+2,nz,nx);

if ~nemflag

ndata=size(B.biots,2);

for data = 1:ndata
    poolidx = B.biots(3,data);
    type = B.biots(4,data);
           if type == -6
               Ctmp = B.biots(5,data);%catch (in t) over next time step
               Ctmp = Ctmp .* 10.^6; % transform catch in mole N
               % distribute overall catch among boxes. Assume that catch is
               % proportional to fish biomass
               C = Ctmp .* bfrac(:,:,poolidx); % nz x nx array (molN)
               C = C./(G.dz.*G.dx.*(G.area/sum(G.dx(1,:)))); %molN.m-3 nz x nx array
               C = C./G.dt; %molN.m-3.s-1 nz x nx array
               
               F = fishingmortality(C, bv(:,:,poolidx)); % nz x nx array
               fmort(poolidx,B.idx.fish,:,:) = permute(F.* bv(:,:,poolidx),...
                   [4 3 1 2]); % molN.m-3.s-1
               
           elseif type == 4
               F=B.biots(5,data);%s-1 scalar
               fmort(poolidx,B.idx.fish,:,:) = permute(F.*bv(:,:,poolidx),...
                   [4 3 1 2]);
    
           end
    
end
end

function F = fishingmortality(c, biomass)
F= -ln(1-(c./biomass));
   