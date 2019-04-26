function [egest, excrete, graze] = egeexc(pred, graze, B, nb, nz, nx)

alleat = pred + graze; % (never overlap) nb+2 x nb+2 x nz x nx
alleat = permute(sum(alleat,1), [3 4 2 1]);  %  nz x nx x nb+2, total eaten by each predator

[egest, excrete] = deal(zeros(nb+2,nb+2,nz,nx));

% Egestion and excretion

egetmp = bsxfun(@times, alleat, reshape([B.gs' 0 0],1,1,[]));
exctmp = bsxfun(@times, alleat, reshape([(1 - B.gs' - B.ge') 0 0],1,1,[]));

egest(:,B.idx.pon,:,:) = permute(egetmp,[3 4 1 2]);
excrete(:,B.idx.nh4,:,:) = permute(exctmp,[3 4 1 2]);
    
% Assume no critters can assimilate Si, all egested immediately to opal

grasi = graze(B.idx.pl,:,:,:) + pred(B.idx.pl,:,:,:); % 1 x nb+2 x nz x nx
grasi = permute(grasi, [3 4 2 1]);                  %  nz x nx x nb+2 
grasi = sum(grasi,3) .* B.RSiN;                      % nz x nx total amount of plSi eaten in each box

egest(B.idx.plsi,B.idx.opal,:,:) = reshape(grasi,1,1,nz,nx);

