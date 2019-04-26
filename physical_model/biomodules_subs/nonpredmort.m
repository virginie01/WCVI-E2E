function mort = nonpredmort(nemflag, bv, P, B, nb, nz, nx)

mort = zeros(nb+2,nb+2,nz,nx);

% N

if nemflag
 
% attention seems that mor0 and Kmor weren't passed from Np to Biovars in setwcenemparams...
mor0=repmat(reshape(B.mor0,1,1,[]),nz,nx,1);
temp=repmat(P.T,1,1,length(B.mor0));
kmor=repmat(reshape(B.Kmor,1,1,[]),nz,nx,1); 
mortn = bsxfun(@times, bv.^2, tempdep(mor0, kmor, temp)); % nz x nx x nbsv
mort(1:nb,B.idx.pon,:,:) = permute(mortn, [3 4 1 2]);

else
    
mortn = bsxfun(@times, bsxfun(@power, bv, reshape(B.m0exp,1,1,[])), reshape(B.m0coef,1,1,[]));%nz x nx x nbsv
mort(1:nb,B.idx.pon,:,:) = permute(mortn,[3 4 1 2]); % N flux (per critter -> PON)

end

% Si

mort(B.idx.plsi,B.idx.opal,:,:) = permute(mortn(:,:,B.idx.pl) .* B.RSiN, [4 3 1 2]);

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);