function nit = nitrif(bv, G, P, B, nb, nz, nx)
%DENIT denitrification for wce based on Martin Schmidt and Anja Eggert, 2012
% "A regional 3D coupled ecosystem model of the Benguela upwelling system"

dz = G.dz;

% Calculate light limitation factor

shadingbio = sum(bv(:,:,[B.idx.ps B.idx.pl]), 3);

nearsurf = [0.5 0.5];    % m, Extend near to surface (since right at surface divides by 0)
zedge = [nearsurf; cumsum(dz)];
pintedge = [shadingbio(1,:).*nearsurf; cumsum(shadingbio).*dz];
kppedge = B.alpha2 .* pintedge./zedge;

kappaP=zeros(size(zedge,1),nx);
for i=1:nx
kappaP(:,i) = interp1(zedge(:,i), kppedge(:,i), -G.z(:,i));
end

kappa = B.alpha1 + kappaP;

% calculate light in each box

I = P.par24 .* exp(-kappa .* -G.z); %nz x nx

% Inititalise nit

nit = zeros(nb+2,nb+2,nz,nx);

tf = strcmp(B.flux{:,1},'nit');
% [source sink] for denitrification 
nitidx = cell2mat(B.flux{tf,[2 3]}); 

for i=1:size(nitidx,1)
src = nitidx(i,1);
snk = nitidx(i,2);

td = tempdep(B.Nit0, B.KNit, P.T);% nz x nx
ld = lightdep(I, B.I0, B.KI); % nz x nx

nit(src,snk,:,:) = permute(td.*ld.*bv(:,:,src), [4 3 1 2]);

end

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);

%------------------------
% light dependence
%------------------------

function ld = lightdep(I, I0, KI)
ld = (1 - max(0, (I-I0)./(KI+I-I0)));

