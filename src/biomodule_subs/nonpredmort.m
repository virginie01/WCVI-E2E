function mort = nonpredmort(nemflag, bv, basum1, basum2, bfrac, P, B, G, nb, nz, nx)

mort = zeros(nb+2,nb+2,nz,nx);

% N

if nemflag
 
% attention seems that mor0 and Kmor weren't passed from Np to Biovars in setwcenemparams...
mor0=repmat(reshape(B.mor0,1,1,[]),nz,nx,1); % nz x nx x nbsv group-dependent
temp=repmat(P.T,1,1,length(B.mor0)); % nz x nx x nbsv location-dependent
kmor=repmat(reshape(B.Kmor,1,1,[]),nz,nx,1); % nz x nx x nbsv group-dependent 
mortn = bsxfun(@times, bv.^2, tempdep(mor0, kmor, temp)); % nz x nx x nbsv
mort(1:nb,B.idx.pon,:,:) = permute(mortn, [3 4 1 2]);

else
%not including temperature effect for now. It's coherent with EwE as it
%is now. COuld try to include temperature effect later on

%previous version
%mortn = bsxfun(@times, bsxfun(@power, bv, reshape(B.m0exp,1,1,[])), reshape(B.m0coef,1,1,[]));%nz x nx x nbsv
%mort(1:nb,B.idx.pon,:,:) = permute(mortn,[3 4 1 2]); % N flux (per critter -> PON)

%new version- to accomodate quadartic mortality-
% use of density in mol N.m-2 because m0coeff is in s-1.(molN.m-2)-1
% if quadratic shape is preferred because m0coeff has to be taken from EwE-
% Using an average concentration for planktonic groups woouldn't work, it
% would be too low-
% total density based on basum1 for nektonic groups and basum 2 for
% planktonic groups- basum2 is preferred for planktonic groups to convert 
% back to concentrations afterwards- in this case total mortality is lower
% than with basum1 but I divide by a lower depth on ULslope so should give
% appropriate mortalities in molN.m-3.s-1 in each box. 

mortn = zeros(nz,nx,nb);

mortnek = bsxfun(@times,bsxfun(@power,basum1,reshape(B.m0exp,1,[])),reshape(B.m0coef,1,[])); % 1 x nb mortality mol N.m-2.s-1
mortnek = mortnek.*G.area;% 1 x nb mortality in mol N.s-1
isnek = reshape(B.isnek,1,1,[]);
mortn(3,1,isnek) = reshape(mortnek(B.isnek')./(G.dz(3,1).*G.dx(3,1).*(G.area./(G.dx(1,1)+G.dx(1,2)))),1,1,[]); % nz x nx x nb molN.m-3.s-1

mortplankton = bsxfun(@times,bsxfun(@power,basum2,reshape(B.m0exp,1,[])),reshape(B.m0coef,1,[])); % 1 x nb mortality mol N.m-2.s-1
mortplankton = mortplankton.*G.area;% 1 x nb mortality in mol N.s-1
mortplankton = bsxfun(@times, bfrac, reshape(mortplankton,1,1,[])); %nz x nx x nb mortality in molN.s-1 in each box

vpb = G.dz.*G.dx.*(G.area./(G.dx(1,1)+G.dx(1,2))); %nz x nx m3
vpb = repmat(vpb,1,1,nb);%nz x nx x nb volume per box m3
vpb(2,2,[1 2 3]) = (990-100).*G.dx(2,2).*(G.area./(G.dx(1,1)+G.dx(1,2)));
vpb(2,2,4) = (990-75).*G.dx(2,2).*(G.area./(G.dx(1,1)+G.dx(1,2)));
vpb(2,2,5) = 50.*G.dx(2,2).*(G.area./(G.dx(1,1)+G.dx(1,2)));
%vpb(3,2,5) = 0; to remove to avoid NaN thereafter

mortplankton = mortplankton./vpb; % nz x nx x nb mortality in molN.m-3.s-1

isplankton = B.isphy | B.iszoo;
isplankton = reshape(isplankton,1,1,[]);

mortn(:,:,isplankton) = mortplankton(:,:,isplankton);% nz x nx x nb mortality molN.m-3.s-1

mort(1:nb,B.idx.pon,:,:) = permute(mortn, [3 4 1 2]);

end

% Si

mort(B.idx.plsi,B.idx.opal,:,:) = permute(mortn(:,:,B.idx.pl) .* B.RSiN, [4 3 1 2]);

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);