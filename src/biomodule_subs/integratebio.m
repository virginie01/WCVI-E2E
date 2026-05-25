function [newbio, dbdt, Splitdbdt, Diag, badthings] = integratebio(fun, nemflag, G, oldbio, P, B, Arch, it, varargin)
% INTEGRATEBIO Integrate biology state variables in WCVIE2E model
%
% [newbio, dbdt, Splitdbdt, Diag, badthings] = integratebio(fun, nemflag, G, oldbio, P, B, Arch, it, solver1, solver2, ...)
%
%   This function integrates the biological state variables in the WCVIE2E
%   model using one or several ODE solvers in sequence. If a solver fails
%   (negative, NaN, or Inf concentrations), the next solver is attempted.
%
% -------------------------------------------------------------------------
% INPUTS:
%
%   fun      : Function handle for the ODEs of biological processes.
%               Must have the signature:
%               [db, Splitdb, Diag] = fun(nemflag, t, oldbio, G, P, B, Arch, it)
%
%   nemflag  : Logical. TRUE if NEMURO-based mode is active.
%
%   G        : Structure with grid and time step parameters:
%               - t  : current model time (s)
%               - dt : time step (s)
%
%   oldbio   : [nz x nx x nbsv] array of current biological concentrations (mol N m^-3)
%
%   P, B, Arch : Structures with model parameters, biological properties, and archive info
%
%   it       : Current time step index
%
%   varargin : List of solver names to try sequentially, e.g.:
%              {'euler', 'ode4', 'ode45', 'implicit'}
%
% -------------------------------------------------------------------------
% OUTPUTS:
%
%   newbio    : [nz x nx x nbsv] biological concentrations after integration
%
%   dbdt      : [nz x nx x nbsv] instantaneous rates of change (mol N m^-3 s^-1)
%
%   Splitdbdt : Structure with dB/dt contributions by process type
%
%   Diag      : Diagnostics returned by ODE solver
%
%   badthings : Logical array [nz x nx x nbsv]. TRUE if integration failed
%               (negative, NaN, or Inf values)
%
% -------------------------------------------------------------------------
% Notes:
%   - If forced biomasses are defined in B.fish, they are applied after each
%     solver step using applyBiomassForcing().
%   - The function automatically switches solvers when instability occurs.
% -------------------------------------------------------------------------
%
% This file was initially derived from the integratebio routine
% developed by Kelly Kearney for the WCE/NEMURO framework and
% subsequently extended and restructured for the WCVI-E2E coastal
% upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Major redevelopment and extensions by Virginie Bornarel (2017–2026)
% include:
%   - adaptation from 1D to 2D biological state integration
%   - integration of fishing mortality and biomass forcing
%   - support for coupled WCVI-E2E/NEMURO configurations
%   - revised solver interfaces and diagnostics
%   - expanded flux and archive handling
%   - updated stability and error-control procedures
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

% === Initialization ===

t = G.t;
dt = G.dt;
[nz,nx,nb] = size(oldbio);
nmethod = length(varargin);

count = 1;

while count <= nmethod
    
    switch varargin{count}
        
        case 'euler'
            
            [dbdt, Splitdbdt, Diag] = feval(fun, nemflag, t, oldbio, G, P, B, Arch, it);
            newbio = oldbio + dbdt.*dt;
            
         if ~nemflag
             
            %force biomass when necessary
            if any(cell2mat(B.fish(2,:))==-1) % if forcing biomasses exist
                 idx = find(cell2mat(B.fish(2,:)) == -1); %find column indices where those forcing biomasses are
                 poolidx = cell2mat(B.fish(1,idx)); % find corresponding critter indices
                 for i = 1:length(poolidx)
                     if poolidx(i)>10 && poolidx(i)<54%nektonic groups
                         newbio(3,1,poolidx(i)) = B.fish{3,idx(i)}.*0.001885.*G.area; % from t.km-2 to total molesN. 
                         newbio(3,1,poolidx(i)) = newbio(3,1,poolidx(i))./(G.dz(3,1).*G.dx(3,1).*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                     else
                         newbio(:,:,poolidx(i)) = B.fish{3,idx(i)};
                     end
                 end    
            end
            
         end
         
            failflag = false;
          
            setappdata(0, 'solvercheck', 1); % debugging
            
            
        case 'ode4'
            
             [newbio, dbdt, Splitdbdt, Diag, failflag] = ode4splitsnonneg(fun,nemflag,[t t+dt],oldbio, G, P, B, Arch, it);
              newbio = newbio(:,:,:,2);
            
            if ~nemflag
            %force biomass when necessary
            if any(cell2mat(B.fish(2,:))==-1) % if forcing biomasses exist
                 idx = find(cell2mat(B.fish(2,:)) == -1); %find column indices where those forcing biomasses are
                 poolidx = cell2mat(B.fish(1,idx)); % find corresponding critter indices
                 for i = 1:length(poolidx)
                     if poolidx(i)>10 && poolidx(i)<54%nektonic groups
                         newbio(3,1,poolidx(i)) = B.fish{3,idx(i)}.*0.001885.*G.area; % from t.km-2 to total molesN. 
                         newbio(3,1,poolidx(i)) = newbio(3,1,poolidx(i))./(G.dz(3,1).*G.dx(3,1).*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                     else
                         newbio(:,:,poolidx(i)) = B.fish{3,idx(i)};
                     end
                 end    
            end
            
            end
            
            setappdata(0, 'solvercheck', 2); % debugging
            
        case 'ode45'
            
            Opt = odeset('nonnegative', 1:numel(oldbio));
            [tout, newbio, dbdt, Splitdbdt] = odewrap(@ode45, fun, nemflag, [t t+dt], oldbio, Opt, G, P, B, Arch, it);  
            failflag = false;
            
            newbio = endonly(newbio);
            
            dbdt = (newbio - oldbio)./dt;
            
            Splitdbdt = Splitdbdt(1);
            
            Diag = [];
            
            if ~nemflag
                
            %force biomass when necessary
            if any(cell2mat(B.fish(2,:))==-1) % if forcing biomasses exist
                 idx = find(cell2mat(B.fish(2,:)) == -1); %find column indices where those forcing biomasses are
                 poolidx = cell2mat(B.fish(1,idx)); % find corresponding critter indices
                 for i = 1:length(poolidx)
                     if poolidx(i)>10 && poolidx(i)<54%nektonic groups
                         newbio(3,1,poolidx(i)) = B.fish{3,idx(i)}.*0.001885.*G.area; % from t.km-2 to total molesN. 
                         newbio(3,1,poolidx(i)) = newbio(3,1,poolidx(i))./(G.dz(3,1).*G.dx(3,1).*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                     else
                         newbio(:,:,poolidx(i)) = B.fish{3,idx(i)};
                     end
                 end    
            end
                        
            end
            
            setappdata(0, 'solvercheck', 3); % debugging
            
        case 'implicit'
            
            [dbdt, Splitdbdt, Diag] = feval(fun, nemflag, t, oldbio, G, P, B, Arch, it);
            
            splt = struct2cell(Splitdbdt);
            
%             nb = size(oldbio,2);
%             for ifd = 1:length(splt)
%                 
%                 fluxin  = squeeze(sum(splt{ifd}, 1))';
%                 fluxout = squeeze(sum(splt{ifd}, 2))';
%                 dbdtpart{ifd} = fluxin(:,1:nb) - fluxout(:,1:nb);
%                 
%             end
            
%             dbdtexplicit = sum(cat(3, dbdtpart{:}), 3);
            
%             lambda = dbdtimplicit./oldbio;
            lambda = dbdt./oldbio; % s-1 nz x nx x nbsv
            lambda(oldbio == 0) = 0;
            
%           newbio = (oldbio + dbdtexplicit.*dt)./(1 - lambda.*dt);
            newbio = oldbio./(1 - lambda.*dt); % nz x nx x nbsv molN.m-3
            
            if ~nemflag 
                
            %force biomass when necessary
            if any(cell2mat(B.fish(2,:))==-1) % if forcing biomasses exist
                 idx = find(cell2mat(B.fish(2,:)) == -1); %find column indices where those forcing biomasses are
                 poolidx = cell2mat(B.fish(1,idx)); % find corresponding critter indices
                 for i = 1:length(poolidx)
                     if poolidx(i)>10 && poolidx(i)<54%nektonic groups
                         newbio(3,1,poolidx(i)) = B.fish{3,idx(i)}.*0.001885.*G.area; % from t.km-2 to total molesN. 
                         newbio(3,1,poolidx(i)) = newbio(3,1,poolidx(i))./(G.dz(3,1).*G.dx(3,1).*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                     else
                         newbio(:,:,poolidx(i)) = B.fish{3,idx(i)};
                     end
                 end    
            end
           
            end
            % TODO note that splits no longer match total for the implicit
            % ones
            
            fail = isnan(newbio) | isinf(newbio) | newbio < 0;
            failflag = any(fail(:));
%             if failflag
%                 disp('testing');
%             end

            setappdata(0, 'solvercheck', 4); % debugging
            
    end
    
    badthings = isnan(newbio) | isinf(newbio) | newbio < 0;
    
%     if t == 249372000
%         disp('Debug stop');
%     end
    
    if failflag || any(badthings(:))
        count = count + 1;
    else
        break;
    end

end

        
    
