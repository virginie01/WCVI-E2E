function nit = nitrif(bv, I, P, B, nb, nz, nx)
%DENIT denitrification for wce based on Martin Schmidt and Anja Eggert, 2012
% "A regional 3D coupled ecosystem model of the Benguela upwelling system"

% Inititalise nit

nit = zeros(nb+2,nb+2,nz,nx);

tf = strcmp(B.flux(:,1),'nit');
% [source sink] for denitrification 
nitidx = cell2mat(B.flux(tf,[2 3])); 

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

%------------------------
% light dependence
%------------------------

function ld = lightdep(I, I0, KI)
ld = (1 - max(0, (I-I0)./(KI+I-I0)));

