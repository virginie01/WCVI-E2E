function [names, nbsv, A] = setstatevars(BioIn)
%SETSTATEVARS State variable setup for mixed_layer wce module
%
% OUTPUTS
%
% nbsv:        scalar. Number of state variables. Corresponds to the number of functional
%              groups in the Ecopath model and the NEMURO variables that don't
%              overlap with Ecopath. Two functional groups (ZL1 and ZL2) are added 
%              if diapause is turned on.
%
% names:       nbsv x 3 cell array. Column 1 = short name; column 2 = long name;
%              column 3 = units. Same order as in Ecopath and NEMURO variables that don't overlap 
%              with Ecopath are placed at the end in same order as in nemnames. Short names are the 
%              same as in type. Long names are the ones in Ecopath (for NEMURO variables that don't 
%              overlap with Ecopath = nemnames(:,2)). 
%              ZL1 and ZL2 are added at the end if diapause is on.
%
% A:           Structure with the following fields:
%              
%              idx: location/index of NEMURO variables. Two groups
%                   "mystery" and "fisheries" are added at the end of the 
%                    list.
%
%              nemidx: column vector indicating index/location of NEMURO variables in
%                      the cell array {names}. The length of the vector depends on the
%                      number of NEMURO variables.
%
%              extrazooidx: column vector indicating index/location of extrazoo 
%                           variables in the cell array {names}. The length of this 
%                           vector depends on the number of extrazoo variables. Optional 
%                           output (i.e. only when BioIn.isnem is turned off)
%
%              nekidx: column vector indicating index/location of nekton variables in
%                      the cell array {names}. The length of this vector depends on
%                      the number of nekton variables.Optional output (i.e. only
%                      when BioIn.isnem is turned off)

%              is[extrazoo]/[nek]/[zoo]/[phy]/[det]: logical vector
%
%              links: nbsv x nbsv matrix. Code used to describe different
%                     trophic interactions (see lines 134-140) and 0 when no trophic
%                     interaction takes place according to the DC matrix in
%                     Ecopath. Note: if diapause is on, ZL1 and ZL2 don't eat and
%                     don't get eaten.
%  
%              fisheries: nbsv x ngear matrix. Contains 1 when species are
%                         harvested by a given gear and 0 otherwise.
%
%              nekdist: Imported from BioIn.nekdist
%
%
%
%-------------------------
% Set state variable names
%-------------------------

nemnames = {...
    'PS',       'PS',          'molN/m^3'
    'PL',       'PL',          'molN/m^3'
    'ZS',       'ZS',          'molN/m^3'
    'ZL',       'ZL',          'molN/m^3'
    'ZP',       'ZP',          'molN/m^3'
    'PON',      'PON',         'molN/m^3'
    'NO3',      'NO3',         'molN/m^3'
    'NH4',      'NH4',         'molN/m^3'
    'DON',      'DON',         'molN/m^3'
    'SiOH4',    'SiOH4',       'molSi/m^3'
    'Opal',     'Opal',        'molSi/m^3'
    'PLsi'      'PLsi',        'molSi/m^3'};

if BioIn.isnem
    names = nemnames;
    [tf, loc] = ismember(lower(nemnames(:,1)), BioIn.types); 
    names(tf,2) = BioIn.EM.name(loc(tf));
    nbsv = size(names,1);
else
    % Figure out which of the Ecopath model groups correspond to these
    % nemuro-derived variables.

    [tf, loc] = ismember(lower(nemnames(:,1)), BioIn.types); 

    nbsv = BioIn.EM.ngroup + sum(~tf);

    names = cell(nbsv,3);
    
    names(loc(tf),:) = nemnames(tf,:);
    names(loc(tf),2) = BioIn.EM.name(loc(tf));
    names(BioIn.EM.ngroup+1:end,:) = nemnames(~tf,:);

    % Name the remaining nekton and zooplankton groups, using N# and Z# for
    % short names and Ecopath names for long names

    isemp = cellfun('isempty', names(:,2));
    names(isemp,2) = BioIn.EM.name(isemp);

    [names{isemp,3}] = deal('mol N m^-3');
    isn = strcmp(BioIn.types, 'n');
    isz = strcmp(BioIn.types, 'z');

    names(isn,1) = cellstr(num2str((1:sum(isn))', 'N%02d'));
    names(isz,1) = cellstr(num2str((1:sum(isz))', 'Z%02d'));
end

%-------------------------
% Identifiers
%-------------------------

A.idx.ps    = find(strcmp(names(:,1), 'PS'));
A.idx.pl    = find(strcmp(names(:,1), 'PL'));
A.idx.no3   = find(strcmp(names(:,1), 'NO3'));
A.idx.nh4   = find(strcmp(names(:,1), 'NH4'));
A.idx.sioh4 = find(strcmp(names(:,1), 'SiOH4'));
A.idx.don   = find(strcmp(names(:,1), 'DON'));
A.idx.pon   = find(strcmp(names(:,1), 'PON'));
A.idx.opal  = find(strcmp(names(:,1), 'Opal'));
A.idx.plsi  = find(strcmp(names(:,1), 'PLsi'));
A.idx.zs    = find(strcmp(names(:,1), 'ZS'));
A.idx.zl    = find(strcmp(names(:,1), 'ZL'));
A.idx.zp    = find(strcmp(names(:,1), 'ZP'));


%-------------------------
% Links
%-------------------------

if ~BioIn.isnem
    
    % Marker for nemuro-derived variables, extra zoo, and nekton

    [~, A.nemidx] = ismember(nemnames(:,1), names(:,1));
    [~, A.extrazooidx]=ismember(names(isz,1),names(:,1));
    [~, A.nekidx]=ismember(names(isn,1),names(:,1));
    
else
    A.nemidx = 1:nbsv;
end
    
    %nolive= or(strcmp(names(:,1), 'N03'), strcmp(names(:,1),'NH4'),...
    %    strcmp(names(:,1),'SiOH4'),strcmp(names(:,1),'DON'),...
    %    strcmp(names(:,1),'PON'),strcmp(names(:,1),'Opal'),...
    %    strcmp(names(:,1),'PLsi'));
    
    %[~,liveidx]=ismember(names(~nolive,1),names(:,1));
    
if ~BioIn.isnem
    A.isextrazoo = regexpfound(names(:,1), 'Z[0-9][0-9]');
    A.isnek = regexpfound(names(:,1), 'N[0-9][0-9]');
    A.iszoo = ismember(names(:,1), {'ZS','ZL','ZP'}) | A.isextrazoo;
    A.isphy = ismember(names(:,1), {'PS','PL'});
    A.isdet = ismember(names(:,1), {'PON'});
else
    A.iszoo = ismember(names(:,1), {'ZS','ZL','ZP'});
    A.isphy = ismember(names(:,1), {'PS','PL'});
    A.isdet = ismember(names(:,1), {'PON'});   
end
    

    % Links
    
if ~BioIn.isnem 
    
        A.links = zeros(nbsv);
        A.links(bsxfun(@and, A.iszoo, A.iszoo')) = 1; % Z-eat-Z
        A.links(bsxfun(@and, A.iszoo, A.isnek')) = 2; % N-eat-Z
        A.links(bsxfun(@and, A.isnek, A.isnek')) = 3; % N-eat-N
        A.links(bsxfun(@and, A.isphy, A.iszoo')) = 4; % Z-eat-P
        A.links(bsxfun(@and, A.isphy, A.isnek')) = 5; % N-eat-P
        A.links(bsxfun(@and, A.isdet, A.iszoo')) = 6; % Z-eat-Det
        A.links(bsxfun(@and, A.isdet, A.isnek')) = 7; % N-eat-Det
        [r,c] = find(table2array(BioIn.EM.dc) == 0);
        idx = sub2ind(size(A.links), r, c);
        A.links(idx) = 0;
    
        A.fisheries = zeros(nbsv,BioIn.EM.ngear);
        [r,c] = find(table2array(BioIn.EM.landing)~=0 | table2array(BioIn.EM.discard)~=0);
        idx = sub2ind(size(A.fisheries), r, c);
        A.fisheries(idx) = 1;
    
else
        A.links = zeros(nbsv);
        A.links(bsxfun(@and, A.iszoo, A.iszoo')) = 1; % Z-eat-Z
        A.links(bsxfun(@and, A.isphy, A.iszoo')) = 4; % Z-eat-P
        A.links(bsxfun(@and, A.isdet, A.iszoo')) = 6; % Z-eat-Det
        
        nemvar = ismember(BioIn.types, {'ps';'pl';'zs';'zl';'zp';'pon'});
        nemdc = table2array(BioIn.EM.dc);
        nemdc = nemdc(nemvar,nemvar);
        
        preypred = ismember(lower(names(:,1)),{'ps';'pl';'zs';'zl';'zp';'pon'});
        tmplinks = A.links(preypred,preypred);
        
        [r,c] = find(nemdc==0);
        idx = sub2ind(size(nemdc),r,c);
        tmplinks(idx) = 0;
        
        A.links(preypred,preypred) = tmplinks;

end       
%-------------------------
% Diapause-related
%-------------------------

A.diapause = BioIn.diapause;
if A.diapause
    names = [
        names
        {...
        'ZL1',       'LargeZooplanktonNoDiapause', 'mol N m^-3'
        'ZL2',       'LargeZooplanktonDiapause',   'mol N m^-3'}];
    A.idx.zl1 = nbsv+1;
    A.idx.zl2 = nbsv+2;
    
    nbsv = nbsv + 2;
 
    linktmp = zeros(nbsv);
    linktmp(1:nbsv-2,1:nbsv-2) = A.links;
    A.links = linktmp;
        
    if ~BioIn.isnem
        A.isnek = [A.isnek; false; false];
        
        linktmp = zeros(nbsv,BioIn.EM.ngear);
        linktmp(1:nbsv-2,1:BioIn.EM.ngear)=A.fisheries;
        A.fisheries = linktmp;
        
    end

end

%-------------------------
% Non-state-variable 
% identifiers
%-------------------------

%if ~BioIn.isnem
A.idx.mys   = nbsv + 1; % "Mystery" box, for source/sink terms coming from or going to nowhere
A.idx.fish  = nbsv + 2; % Fisheries

%isnek = [A.isnek;false;false];
%iszoo = [A.iszoo;false;false];

%isimp = false([nbsv+2 1]);
%isimp(A.idx.mys) = true;

%linktmp = zeros(nbsv+2);
%linktmp(1:nbsv, 1:nbsv) = A.links;
%A.links = linktmp;
%A.links(bsxfun(@and,isimp,iszoo')) = 8; %Z-eat-import
%A.links(bsxfun(@and,isimp,isnek')) = 9; %N-eat-import

%c = BioIn.EM.groupdata.import == 0;
%c = c';
%c = [c,true([1 nbsv+2-length(c)])];

%A.links(isimp,c) = 0;

%end



    
