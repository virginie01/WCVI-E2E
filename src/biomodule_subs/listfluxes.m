function list = listfluxes(type, ismixed, Idx, diapause, varargin)
%LISTFLUXES Returns predator-prey and biogeochemical flux definitions
%
%   list = listfluxes('nemuro', ismixed, Idx, diapause)
%   list = listfluxes('WCVI-E2E', ismixed, Idx, diapause, links, fisheries)
%
%   DESCRIPTION:
%   Constructs a list of ecological and biogeochemical fluxes for the
%   WCVI-E2E end-to-end model. Depending on the 'type' specified, it
%   returns fluxes for either the NEMURO-style ('nemuro') or
%   WCVI-E2E ('WCVI-E2E') model configuration.
%
% INPUT VARIABLES:
%
% type    : String. Either:
%               - 'nemuro' → NEMURO model formulation
%               - 'WCVI-E2E' → WCVI-E2E coupled model formulation
%
% ismixed : Logical vector indicating which state variables are
%                 subject to vertical or horizontal mixing.
%
% Idx     : Structure containing indices of NEMURO state variables (e.g.,
%           Idx.ps, Idx.pl, Idx.zs, Idx.no3, etc.).
%
% diapause: Logical scalar. If true, includes vertical fluxes for
%                 diapausing zooplankton (e.g., ZL1).
%
% links   : [Only for 'wce'] ng x ng predator-prey interaction matrix. 
%           Each element indicates type of interaction:
%           1 = zooplankton eat zooplankton
%           2 = nekton eat zooplankton
%           3 = nekton eat nekton
%           4 = zooplankton eat phytoplankton
%           5 = Nekton eat phytoplankton
%           6 = Zooplankton eat detritus
%           7 = Nekton eat detritus
%
% fisheries   : [Only for 'wce'] ng x ngear matrix.
%               1 if the species is harvested by a gear, 0 otherwise.
%
% OUTPUT VARIABLE:
%
%   list        : n x 3 cell array, where each row describes one flux:
%                   Column 1 → flux type (string)
%                   Column 2 → source index
%                   Column 3 → sink index
%
% NOTES:
%   - Assumes all “living” groups are included in at least one predator–
%     prey interaction. If not, mortality fluxes may need updating.
%   - Physical flux placeholders (V, H, CS, etc.) are added for completeness.
%
% This file was initially derived from the original listfluxes
% routine developed by Kelly Kearney for the WCE/NEMURO framework
% and subsequently extended and restructured for the WCVI-E2E
% coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2011 Kelly Kearney
%
% Major modifications and extensions by Virginie Bornarel (2017–2026)
% include:
%   - expanded trophic interaction categories
%   - fisheries and physical transport flux bookkeeping
%   - support for diapause and mixed-variable handling
%   - separate predator inflow/outflow accounting
%   - generalized flux-definition architecture for WCVI-E2E
%   - revised documentation and ecosystem process organization
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% -------------------------------------------------------------------------
% Parse optional inputs
% -------------------------------------------------------------------------
if ~isempty(varargin)
links = varargin{1};
fisheries = varargin{2};
end

% -------------------------------------------------------------------------
% Define flux lists by model type
% -------------------------------------------------------------------------
switch type

    %======================================================================
    % NEMURO MODEL ("nemuro")
    %======================================================================
    case 'nemurokak'

    % -------------------------
    % Primary production
    % ------------------------- 
        
        A.gpp = [...
            Idx.no3     Idx.ps
            Idx.nh4     Idx.ps
            Idx.no3     Idx.pl
            Idx.nh4     Idx.pl
            Idx.sioh4   Idx.plsi];

     % -------------------------
     % Exudation and respiration
     % -------------------------

        A.exx = [...
            Idx.ps      Idx.don
            Idx.pl      Idx.don
            Idx.plsi    Idx.sioh4];

        A.res = [...
            Idx.ps      Idx.no3     
            Idx.ps      Idx.nh4     
            Idx.pl      Idx.no3     
            Idx.pl      Idx.nh4     
            Idx.plsi    Idx.sioh4];  

     % -------------------------
     % Grazing and mortality
     % -------------------------
           
        A.gra = [...
            Idx.ps      Idx.zs
            Idx.ps      Idx.zl
            Idx.pl      Idx.zl
            Idx.pl      Idx.zp
            Idx.zs      Idx.zl
            Idx.zs      Idx.zp
            Idx.zl      Idx.zp];
        
        A.pre_out = [];

        A.mor = [...
            Idx.ps      Idx.pon
            Idx.pl      Idx.pon
            Idx.zs      Idx.pon
            Idx.zl      Idx.pon
            Idx.zp      Idx.pon
            Idx.plsi    Idx.opal];

     % -------------------------
     % Excretion, egestion, remineralization
     % -------------------------
        
        A.exc = [...
            Idx.zs      Idx.nh4
            Idx.zl      Idx.nh4
            Idx.zp      Idx.nh4];
        
        A.ege = [...
            Idx.zs      Idx.pon
            Idx.zl      Idx.pon
            Idx.zp      Idx.pon
            Idx.plsi    Idx.opal];
        


        A.fish = [];

        
        A.amm = [...
            Idx.don     Idx.nh4
            Idx.pon     Idx.nh4
            Idx.pon     Idx.don
            Idx.opal    Idx.sioh4];
        
        A.den = [...
            Idx.don     Idx.nh4
            Idx.pon     Idx.nh4
            Idx.pon     Idx.don
           Idx.opal    Idx.sioh4];
        
        A.nit = [...
            Idx.nh4     Idx.no3];
             
     % -------------------------
     % Physical processes
     % ------------------------- 
        
        mixid = find(ismixed);
        
        A.V = [...
            mixid     mixid];
        
        A.H = [...
            mixid     mixid];
        
        A.X = [...
            mixid     mixid];
        
        A.CS = [...
            mixid     mixid];
        
        A.P = [...
            mixid     mixid];
        
        A.R = [...
            mixid     mixid];

        A.VICC = [...
            mixid     mixid];
        
        A.DC = [...
            mixid     mixid];

        A.SBC = [...
            mixid     mixid];
        
        A.CU = [...
            mixid     mixid];
        
     % -------------------------
     % Vertical sinking or migration
     % -------------------------
        if diapause
        A.vflx = [...
            Idx.pon  Idx.pon
            Idx.opal Idx.opal
            Idx.zl  Idx.zl];
        else
        A.vflx = [...
            Idx.pon  Idx.pon
            Idx.opal Idx.opal];
        end 


    %======================================================================
    % WCVI-E2E MODEL
    %======================================================================    
    case 'wce'

        % -------------------------
        % Primary production
        % -------------------------
        
        A.gpp = [...
            Idx.no3     Idx.ps
            Idx.nh4     Idx.ps
            Idx.no3     Idx.pl
            Idx.nh4     Idx.pl
            Idx.sioh4   Idx.plsi];
        
        A.npp = [...
            Idx.mys     Idx.ps
            Idx.mys     Idx.pl];

        % -------------------------
        % Exudation and respiration
        % -------------------------

        A.exx = [...
            Idx.ps      Idx.don
            Idx.pl      Idx.don
            Idx.plsi    Idx.sioh4];

        A.res = [...
            Idx.ps      Idx.no3     
            Idx.ps      Idx.nh4     
            Idx.pl      Idx.no3     
            Idx.pl      Idx.nh4     
            Idx.plsi    Idx.sioh4];   

        % -------------------------
        % Predator–prey relationships
        % -------------------------
        
        [pry1, prd1] = find(links == 2 | links == 3 | links == 5 | links ==7); % nektonic
        [pry2, prd2] = find(links == 1 | links == 4 | links ==6); % planktonic
        
        A.gra = [...
            pry2        prd2];
        
        A.pre_in = [...
            pry1        prd1];
        
        A.pre_out = [...
            pry1        prd1];

        % -------------------------
        % Excretion and egestion
        % -------------------------
        
        unqprd = unique([prd1; prd2]);
        
        A.exc_in = [...
            unqprd      Idx.nh4.*ones(size(unqprd))];% all predator groups excret towards nh4

        A.exc_out = [...
            unqprd      Idx.nh4.*ones(size(unqprd))];% all predator groups excret towards nh4
 
        
        A.ege_in = [...
            unqprd      Idx.pon.*ones(size(unqprd)) % egestion from all predator groups towards PON
            Idx.plsi    Idx.opal];

        A.ege_out = [...
            unqprd      Idx.pon.*ones(size(unqprd)) % egestion from all predator groups towards PON
            Idx.plsi    Idx.opal];

        % -------------------------
        % Mortality and fisheries
        % ------------------------- 
        
        alllive = unique([pry1; pry2; prd1; prd2]); 
        alllive = alllive(alllive~=Idx.pon & alllive ~= Idx.opal);% Assumes all live guys except det
        
        A.mor = [...
            alllive     Idx.pon.*ones(size(alllive))
            Idx.plsi    Idx.opal];
        
        fishedsp = find(~all(fisheries == 0,2));
        
        A.fish = [...
            fishedsp    Idx.fish.*ones(size(fishedsp))]; 

        % -------------------------
        % Remineralization and nutrient cycling
        % -------------------------
        
        A.amm = [...
            Idx.don     Idx.nh4
            Idx.pon     Idx.nh4
            Idx.pon     Idx.don
            Idx.opal    Idx.sioh4];
        
        A.den = [...
            Idx.don     Idx.nh4
            Idx.pon     Idx.nh4
           %Idx.no3     Idx.mys  %no3 and pon produce nh4 in 1 reaction
            Idx.pon     Idx.don
           Idx.opal    Idx.sioh4];
        
        A.nit = [...
            Idx.nh4     Idx.no3];
                
        % -------------------------
        % Physical processes
        % -------------------------
        
        mixid = find(ismixed);
        
        A.V = [...
            mixid     mixid];
        
        A.H = [...
            mixid     mixid];
        
        A.X = [...
            mixid     mixid];
        
        A.CS = [...
            mixid     mixid];
        
        A.P = [...
            mixid     mixid];
        
        A.R = [...
            mixid     mixid];

        A.VICC = [...
            mixid     mixid];
        
        A.DC = [...
            mixid     mixid];

        A.SBC = [...
            mixid     mixid];
        
        A.CU = [...
            mixid     mixid];
        
        % -------------------------
        % Vertical sinking or migration
        % -------------------------
        if diapause
        A.vflx = [...
            Idx.pon  Idx.pon
            Idx.opal Idx.opal
            Idx.zl  Idx.zl];
        else
        A.vflx = [...
            Idx.pon  Idx.pon
            Idx.opal Idx.opal];
        end            

        
    otherwise
        error('Unrecognized option for biological module');
        
        
end
        
% -------------------------------------------------------------------------
% Convert structure A into standardized list
% -------------------------------------------------------------------------

list = cell(0,3);

fld = fieldnames(A);
for ifld = 1:length(fld)
    nlnk = size(A.(fld{ifld}),1);
    
    list = [list; repmat(fld(ifld), nlnk, 1), num2cell(A.(fld{ifld}))];
end
            
        
