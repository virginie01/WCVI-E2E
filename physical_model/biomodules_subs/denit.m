function den = denit(bv, P, B, nb, nz, nx)
%DENIT denitrification for wce based on Martin Schmidt and Anja Eggert, 2012
% "A regional 3D coupled ecosystem model of the Benguela upwelling system"

den = zeros(nb+2,nb+2,nz,nx);

% delimiter fX
fo2 = fdecomp(B.alphao2,B.o2); % nz x nx
fH2S = fdecomp(B.alphaH2S, B.H2S); % nz x nx
fNo3 = fdecomp(B.alphaNo3, bv(:,:,B.idx.no3)); %nz x nx
fNH3 = decomp(B.alphaNH4, bv(:,:,B.idx.nh4)); %nz x nx

tf = strcmp(B.flux{:,1},'den');
% [source sink] for denitrification 
denidx = cell2mat(B.flux{tf,[2 3]}); 

for i=1:size(denidx,1)
src = denidx(i,1);
snk = denidx(i,2);
if all([src snk] ~= [B.idx.no3 B.idx.mys])
td = tempdep(B.vdec(src,snk), B.Kdec(src,snk), P.T);% nz x nx
den(src,snk,:,:) = permute(td.*(1-fo2).*fNo3.*(1-fNH3.*(1-fH2S)), [4 3 1 2]);
elseif all([src snk] == [B.idx.no3 B.idx.mys])
td = tempdep(B.vdec(B.idx.pon,B.idx.nh4), B.Kdec(B.idx.pon,B.idx.nh4), P.T);% nz x nx
den(src,snk,:,:) = permute(-B.RNo3PON.*td.*(1-fo2).*fNo3.*(1-fNH3.*(1-fH2S)), [4 3 1 2]);
end
end


%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);

%-------------------------------------
% delimiter fX from 
% Martin Schmidt and Anja Eggert, 2012
% "A regional 3D coupled ecosystem model 
% of the Benguela upwelling system"
%-------------------------------------

function fX = fdecomp (alphax, X)
fX = (1- exp(-2.*alphax.*X))./(1+ exp(-2.*alphax.*X));

