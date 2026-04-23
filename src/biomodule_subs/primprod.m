function [gpp, exc, resp, dz, psmax, Lfc, no3lim, nh4lim, silim, I, kappa, kappaP] = primprod(bv, G, P, B, nz, nx, nb)
%PRIMPROD Primary production fluxes for wce module

% OUTPUTS
% 
% gpp:    nbsv x nbsv x nz x nx array. Gross primary production fluxes (source, sink, z domain, x domain)
%         at current time step
%
% exc:    nbsv x nbsv x nz x nx. Extracellular excretion fluxes (source, sink, z domain, x domain)
%         at current time step
%
% resp:   nbsv x nbsv x nz x nx. Phytoplankton respiration fluxes (source, sink, z domain, x domain)
%         at current time step
%
% dz:     nz x nx array. Depth intervals at current time step
%
% psmax:  nz x nx x nbsv array. Maximum possible production rate at current temperature
%
% Lfc:    nz x nx x nbsv array. Light limitation factor at current time step
%
% no3lim: nz x nx x nbsv array. N03 Limitation at current time step
%
% nh4lim: nz x nx x nbsv array. NH4 Limitation at current time step
%
% silim:  nz x nx x nbsv array. SiOH4 Limitation at current time step
%
% I:      nz x nx x nbsv array. Light intensity at current time step which is location-dependent 
%         (i.e. same values along the 3rd dim)
%
% kappa:  nz x nx array. 
%
% kappaP: nz x nx array.


dz = G.dz;

% Calculate light limitation factor

shadingbio = sum(bv(:,:,[B.idx.ps B.idx.pl]), 3); %nz x nx array total phytoplankton concentration molN/m3

%% Kelly's version of kappaP

nearsurf = [0.5 0.5];    % m, Extend near to surface (since right at surface divides by 0)
zedge = [nearsurf; cumsum(dz)]; %(nz+1) x nx array zedge
pintedge = [shadingbio(1,:).*nearsurf; cumsum(shadingbio).*dz]; %(nz+1) x nx array molN.m-2
kppedge = B.alpha2 .* pintedge./zedge; %(nz+1) x nx array m-1

kappaP=zeros(nz,nx);
for i=1:nx
kappaP(:,i) = interp1(zedge(:,i), kppedge(:,i), G.z(:,i)); %nz x nx array
end

%% My version of kappaP

%a = shadingbio.*dz;%nz x nx molN/m2 in each box

%b = [a(1,:)./2; ...
%    a(1,:)+(a(2,:)./2);
%    a(1,:)+a(2,:)+(a(3,:)./2)];%nz x nx molN/m2 at the center of each box

%c = b./G.z;%nz x nx molN/m3 cumul at the center of each box

%kappaP = B.alpha2 .* c; %nz x nx array m-1

kappa = B.alpha1 + kappaP; %nz x nx array m-1

I = zeros(nz,nx);
for i = 1:nx
I(:,i) = P.par24(i) .* exp(-kappa(:,i) .* G.z(:,i)); % nz x nx array W.m-2
end

if B.usesteele
    I = repmat(I,1,1,length(B.Iopt));               %nz * nx * nbsv I is location-dependent 
    Iopt = repmat(ones(nz,nx),1,1,length(B.Iopt));
    Iopt = Iopt.*reshape(B.Iopt,1,1,[]);            %nz x nx x nbsv Iopt is group-dependent
    Lfc = I./Iopt .* exp(1 - I./Iopt);              %nz x nx x nbsv
else
    I = repmat(I,1,1,length(B.alpha));              %nz * nx * nbsv I is location-dependent 
    alpha = repmat(ones(nz,nx),1,1,length(B.alpha));
    alpha = alpha.*reshape(B.alpha,1,1,[]);         %nz x nx x nbsv alpha is group-dependent 
    Vmax = repmat(ones(nz,nx),1,1,length(B.Vmax));
    Vmax = Vmax.*reshape(B.Vmax,1,1,[]);            %nz x nx x nbsv Vmax is group dependent
    Lfc = 1 - exp(-alpha .* I./Vmax);               %nz x nx x nbsv
end

% Calculate gross primary production

[no3lim, nh4lim, silim, fratio, psmax] = deal(zeros(nz, nx, nb));
[gpp, exc, resp] = deal(zeros(nb+2, nb+2, nz, nx));

for ib = [B.idx.ps B.idx.pl]
    
    % Macronutrient limitation
    
    no3lim(:,:,ib) = mich(bv(:,:,B.idx.no3), B.Kno3(ib)) .* exp(-B.pusai(ib) .* bv(:,:,B.idx.nh4));
    nh4lim(:,:,ib) = mich(bv(:,:,B.idx.nh4), B.Knh4(ib));
    silim(:,:,ib)  = mich(bv(:,:,B.idx.sioh4), B.Ksi(ib)); % Always 1 if Ksi == 0, i.e. not si-dependant
    
    % F-ratio (fraction new production)   
    
    fratio(:,:,ib) = no3lim(:,:,ib)./(no3lim(:,:,ib) + nh4lim(:,:,ib));% nz x nx array
      
    % Overal nutrient limitation
    
    % nutlim = min([(no3lim(:,ib) + nh4lim(:,ib)) silim(:,ib) felim(:,ib)], [], 2);
    if ib == B.idx.ps
        nutlim = no3lim(:,:,ib) + nh4lim(:,:,ib); %nz x nx array
    else   
        nutlim = min(cat(3,(no3lim(:,:,ib) + nh4lim(:,:,ib)), silim(:,:,ib)), [], 3); %nz x nx array
    end
    
    % Maximum possible production rate at current temps
    
    psmax(:,:,ib) = tempdep(B.Vmax(ib), B.Kgpp(ib), P.T); %nz x nx array
    
    % Nitrogen uptake
    
    totnuptake = psmax(:,:,ib) .* Lfc(:,:,ib) .* nutlim .* bv(:,:,ib);%nz x nx array
    
    no3uptake = totnuptake .* fratio(:,:,ib); %nz x nx array
    nh4uptake = totnuptake .* (1 - fratio(:,:,ib)); %nz x nx array
    
    % Uptake of SiOH4
    
    if ib == B.idx.pl
        siuptake = (no3uptake + nh4uptake) .* B.RSiN; %nz x nx array
    else
        siuptake = zeros(nz,nx); %
    end
    
    % Gross primary production fluxes
    
    gpp(B.idx.no3,ib,:,:) = no3uptake;
    gpp(B.idx.nh4,ib,:,:) = nh4uptake;
    if ib == B.idx.pl
        gpp(B.idx.sioh4,B.idx.plsi,:,:) = siuptake;
    end

    % Extracellular excretion (N) is proportional to uptake
    
    excn = B.gamma(ib) .* (no3uptake + nh4uptake);
    exc(ib,B.idx.don,:,:) = excn;
        
    % Phytoplankton respiration (N)
    
    respn = tempdep(B.res0(ib), B.Kres(ib), P.T) .* bv(:,:,ib); %nz x nx array
    
    resp(ib,B.idx.no3,:,:) = respn .* fratio(:,:,ib);
    resp(ib,B.idx.nh4,:,:) = respn .* (1 - fratio(:,:,ib));
    
    % Excretion and respiration of silica proportional to N
    
    if ib == B.idx.pl
        excsi = excn * B.RSiN;
        exc(B.idx.plsi,B.idx.sioh4,:,:) = excsi;
        
        respsi = respn * B.RSiN;
        resp(B.idx.plsi, B.idx.sioh4,:,:) = respsi;
    end
    
    % Excretion and respiration of iron proportional to N
    % TODO: with quota model, should Fe be excreted and respired?
    
%     excfe = excn .* (A.RCN/A.RCFe);
%     exc(nb+1,A.idx.fe,:) = permute(excfe, [2 3 1]) + exc(nb+1,A.idx.fe,:);
%     
%     respfe = respn .* (A.RCN/A.RCFe);
%     resp(nb+1,A.idx.fe,:) = permute(respfe, [2 3 1]) + resp(nb+1,A.idx.fe,:);
    
end

%------------------------
% Michaelis-Menten uptake
% limitation
%------------------------

function lim = mich(x, kx)
lim = x./(x + kx);
lim(isnan(lim)) = 0; % if x == 0 and kx == 0

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);
