function list = listfluxes(type, ismixed, Idx, diapause, varargin)
%LISTFLUXES Return indices of fluxes in nemurokak or wce model
%
% list = listfluxes('nemurokak', Idx)
% list = listfluxes('wce', Idx, links)
%
% Input variables:
%
%   Idx:    1 x 1 structure with fields (ps, pl, zs, zl, etc) indicating
%           the index of the biological state variables corresponding to
%           each nemuro-derived variable (see wce and nemurokak setup)
%
%   links:  ng x ng array indicating type of predator/prey interaction
%           1 = zooplankton eat zooplankton
%           2 = nekton eat zooplankton
%           3 = nekton eat nekton
%           4 = zooplankton eat phytoplankton
%           5 = Nekton eat phytoplankton
%           6 = Zooplankton eat detritus
%           7 = Nekton eat detritus
%           Note: For mortality flux indices, this function assumes that
%           all living critters are involved in at least one predator-prey
%           interaction.  If not true, I'll need to update this.
%
% Output variables:
%
%   list:   n x 3 cell array, where column 1 indicates the type of flux,
%           column 2 source group index, and column 3 the sink group index

% Copyright 2011 Kelly Kearney

if ~isempty(varargin)
links = varargin{1};
fisheries = varargin{2};
end

switch type
    
    case 'nemurokak'
        
        A.gpp = [...
            Idx.no3     Idx.ps
            Idx.nh4     Idx.ps
            Idx.no3     Idx.pl
            Idx.nh4     Idx.pl
            Idx.sioh4   Idx.plsi];

        A.exx = [...
            Idx.ps      Idx.don
            Idx.pl      Idx.don
            Idx.plsi    Idx.sioh4];
%             Idx.mys     Idx.fe];

        A.res = [...
            Idx.ps      Idx.no3     
            Idx.ps      Idx.nh4     
            Idx.pl      Idx.no3     
            Idx.pl      Idx.nh4     
            Idx.plsi    Idx.sioh4];   
%             Idx.mys     Idx.fe];
        
        A.gra = [...
            Idx.ps      Idx.zs
            Idx.ps      Idx.zl
            Idx.pl      Idx.zl
            Idx.pl      Idx.zp
            Idx.zs      Idx.zl
            Idx.zs      Idx.zp
            Idx.zl      Idx.zp];
        
        A.pre_out = [];
        
        A.exc = [...
            Idx.zs      Idx.nh4
            Idx.zl      Idx.nh4
            Idx.zp      Idx.nh4];
%           Idx.mys     Idx.fe];
        
        A.ege = [...
            Idx.zs      Idx.pon
            Idx.zl      Idx.pon
            Idx.zp      Idx.pon
            Idx.plsi    Idx.opal];
%             Idx.mys     Idx.fe];
        
        A.mor = [...
            Idx.ps      Idx.pon
            Idx.pl      Idx.pon
            Idx.zs      Idx.pon
            Idx.zl      Idx.pon
            Idx.zp      Idx.pon
            Idx.plsi    Idx.opal];
%             Idx.mys     Idx.fe];

        A.fish = [];

        
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
             
        % terms related to physical fluxes
        
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
        
        %terms related to vertical migration/sinking
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
        
    case 'wce'
        
        A.gpp = [...
            Idx.no3     Idx.ps
            Idx.nh4     Idx.ps
            Idx.no3     Idx.pl
            Idx.nh4     Idx.pl
            Idx.sioh4   Idx.plsi];
        
        A.npp = [...
            Idx.mys     Idx.ps
            Idx.mys     Idx.pl];

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
        
        [pry1, prd1] = find(links == 2 | links == 3 | links == 5 | links ==7); % nektonic
        [pry2, prd2] = find(links == 1 | links == 4 | links ==6); % planktonic
        
        A.gra = [...
            pry2        prd2];
        
        A.pre_in = [...
            pry1        prd1];
        
        A.pre_out = [...
            pry1        prd1];
        
        
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
        
        alllive = unique([pry1; pry2; prd1; prd2]); 
        alllive = alllive(alllive~=Idx.pon & alllive ~= Idx.opal);% Assumes all live guys except det
        
        A.mor = [...
            alllive     Idx.pon.*ones(size(alllive))
            Idx.plsi    Idx.opal];
        
        fishedsp = find(~all(fisheries == 0,2));
        
        A.fish = [...
            fishedsp    Idx.fish.*ones(size(fishedsp))];   
        
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
                
        % terms related to physical fluxes
        
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
        
        %terms related to vertical migration/sinking
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
        
% Reformat into list

list = cell(0,3);

fld = fieldnames(A);
for ifld = 1:length(fld)
    nlnk = size(A.(fld{ifld}),1);
    
    list = [list; repmat(fld(ifld), nlnk, 1), num2cell(A.(fld{ifld}))];
end
            
        
