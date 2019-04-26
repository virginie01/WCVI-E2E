function amm = ammonif(bv, P, B, nb, nz, nx)
%DECOMPREMIN Decomposition, remineralization, etc., for wce

% NEMURO-derived decomposition/remineralization

amm = zeros(nb+2,nb+2,nz,nx);
for isrc = 1:nb
    for isnk = 1:nb
        if B.vdec(isrc,isnk) > 0
            td = tempdep(B.vdec(isrc,isnk), B.Kdec(isrc,isnk), P.T); % nz x nx
            fo2=fdecomp(B.alphao2,B.o2); % nz x nx
            amm(isrc,isnk,:,:) = permute(td .* fo2 .* bv(:,:,isrc), [4 3 1 2]);          
        end
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
