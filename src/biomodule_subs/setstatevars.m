function [names, nbsv, A] = setstatevars(BioIn)
% SETSTATEVARS State variable setup for WCVI-E2E biological module
%
% This function defines and organizes the state variables used in the
% WCVI-E2E biological model (Ecopath-NEMURO hybrid). It returns their
% names, units, count, and indexing information needed for simulation
% of trophic interactions and physical-biogeochemical processes.
%
% INPUT:
% BioIn : Structure holding biological configuration data, including:
% - BioIn.isnem: logical, true if using NEMURO-only mode
% - BioIn.EM: Ecopath model info
% - BioIn.types: group types (e.g. 'z', 'n', etc.)
% - BioIn.diapause: logical, include diapause groups if true
%
% OUTPUTS:
%
%   nbsv       Scalar. Total number of state variables. Includes:
%              - Functional groups defined in the Ecopath model
%              - NEMURO variables not overlapping with Ecopath
%              - (Optional) Two additional groups (ZL1 and ZL2) if diapause is enabled
%
%   names      Cell array of size [nbsv x 3]. Describes each state variable:
%              - Column 1: Short name (e.g., 'ZS', 'NO3'), matching BioIn.types
%              - Column 2: Long name (from Ecopath or nemnames)
%              - Column 3: Units (e.g., 'mol N m^-3')
%              The order is: Ecopath groups first, followed by unmatched NEMURO variables,
%              and ZL1/ZL2 at the end if diapause is active.
%
% A:           Struct containing metadata and mappings related to state variables:
%              
%       A.idx: Index of specific NEMURO or derived variables: e.g.,
%              A.idx.zs, A.idx.no3, A.idx.pon, etc...
%
%       A.nemidx: Vector. Indices of NEMURO variables in the `names` array.
%
%       A.extrazooidx: Vector. Indices of extra zooplankton variables in `names`.
%                       Present only if BioIn.isnem is false.
%
%       A.nekidx: Vector. Indices of nekton variables in `names`.
%                       Present only if BioIn.isnem is false.
%
%       A.is[extrazoo]/[nek]/[zoo]/[phy]/[det]: logical vector
%
%       A.links: [nbsv x nbsv] matrix. Describes trophic interactions:
%                - Value depends on consumer-resource interaction type
%                - Value = 0 if no interaction (based on Ecopath DC matrix)
%                - If diapause is active, ZL1 and ZL2 have no interactions
%  
%       A.fisheries: [nbsv x ngear] matrix. Logical matrix indicating if a 
%                    group is harvested by a specific gear (1 = yes, 0 = no)
%
%
%       A.diapause: Logical. Indicates if diapause mode is active.
%
%       A.idx.mys:  Index of the synthetic "mystery" group (for unaccounted source/sink flows).
%
%       A.idx.fish: Index of the synthetic "fisheries" group (aggregated catch group).
%
%       A.dzmod_fn: The dz modifier function handle (or [] if none). Passed through from BioIn.
%
% This file was derived from the original setstatevars routine
% developed by Kelly Kearney for the WCE/NEMURO framework and
% substantially modified for the WCVI-E2E coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Major modifications and extensions by Virginie Bornarel (2017–2026) include:
%   - revised state-variable ordering for WCVI-E2E/NEMURO configurations
%   - expanded trophic interaction categories, including detrital feeding
%   - addition of fisheries metadata and harvested-group bookkeeping
%   - support for extra zooplankton and nekton index tracking
%   - revised NEMURO-only interaction filtering using Ecopath diet data
%   - updated diapause-group handling and model metadata
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

%------------------------------------------------
% Define standard NEMURO variable names and units
%------------------------------------------------


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

%-------------------------------------
% Initialize variable names and count
%-------------------------------------


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

%--------------------------
% Assign key group indices
%--------------------------

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


%-------------------------------------------
% Determine NEMURO, zoo, nekton identifiers
%-------------------------------------------

if ~BioIn.isnem
    
    % Marker for nemuro-derived variables, extra zoo, and nekton

    [~, A.nemidx] = ismember(nemnames(:,1), names(:,1));
    [~, A.extrazooidx]=ismember(names(isz,1),names(:,1));
    [~, A.nekidx]=ismember(names(isn,1),names(:,1));
    
else
    A.nemidx = 1:nbsv;
end
    
    
    
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
    

%----------------------------------        
% Build trophic interaction matrix
%----------------------------------  
    
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



    
