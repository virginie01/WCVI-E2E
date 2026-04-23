function [newbio, dbdt, Splitdbdt, Diag, badthings] = integratebio(fun, nemflag, G, oldbio, P, B, Arch, it, varargin)
% INTEGRATEBIO Integrate biology in wce biological module
%
% [newbio, dbdt, Splitdbdt, Diag, badthings] = integratebio(fun, G, oldbio, P, B, solver1, solver2, ...)
%
% This function is a wrapper for different ODE solvers, used to integrate biological variables in 
% the wce module.  It allows one to try multiple different solvers, which can sometimes be useful 
% when values get low or start changing quickly (such that they are difficult to solve with
% the set time step). 
%
% Input variables:
%
%   fun:     function handle of ODEs for biology.  Must be of the form
%            [db, Splitdb, Diag] = fun(t, b, ParamsG, ParamsP, ParamsB).
%
%   oldbio:  nz x nx x nbsv array of biological state variables
%
%   G, P, B: structure of additional parameters for ODEs (parameters are for current time step when 
%            applicable)
%
%   solver: ODE solvers to use, in order.  If any tracers go negative, become NaNs, or become 
%           infinite, the next solver will be tried, until the last solver is reached.  Can be 
%           'euler', 'ode4', 'ode45', or 'implicit'.  *NOTE* Implicit doesn't really work,
%           doesn't conserve mass
%
%
% Output variables:
%
%   newbio: biological state variable values at time t+dt
%
%   dbdt:   dB/dt for each state variable.
%
%   Split:  structure with contribution toward dB/dt from each flux type in
%           the ODE function.  For euler and ode4, these will add to dbdt;
%           for ode45 they represent the splits when the ODE function is
%           evaluated at time t (since getting weights throughout the
%           variable steps is not presently possible)
%
%   Diag:   additional diagnostic variables returned by the ODE function.
%           For euler and ode4, the diagnostics are those associated with
%           the beginning of the time step.  Extra diagnostics cannot be
%           returned from the ode45 solver at this time.
%
%   bad:    nz x nx x nbsv logical array, true if the final solver attempted
%           still failed to integrate without hitting a negative, NaN, or
%           Inf value.

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

        
    
