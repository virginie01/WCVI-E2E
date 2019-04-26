function [newbio, dbdt, Splitdbdt, Diag, badthings] = integratebio(fun, G, oldbio, P, B, varargin)
%INTEGRATEBIO Integrate biology in mixed_layer module
%
% [newbio, dbdt, Split, Diag, bad] = integratebio(fun, t, dt, oldbio, ...
%                                    Param, solver1, solver2, ...)
%
% This function is a wrapper for different ODE solvers, used to integrate
% biological variables in a few mixed_layer modules.  It allows one to try
% multiple different solvers, which can sometimes be useful when values get
% low or start changing quickly (such that they are difficult to solve with
% the set time step). 
%
% Input variables:
%
%   fun:    function handle of ODEs for biology.  Must be of the form
%           [db,Splitdb,Diag] = fun(t,b,P).
%
%   t:      current time of integration
%
%   dt:     time step
%
%   oldbio: nz x nx x nbsv array of biological state variables
%
%   P, B:   structure of additional parameters for ODEs
%
%   solver: ODE solvers to use, in order.  If any tracers go negative,
%           become NaNs, or become infinite, the next solver will be tried,
%           until the last solver is reached.  Can be 'euler', 'ode4', or
%           'ode45', or 'implicit'.  *NOTE* Implicit doesn't really work,
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
%   bad:    nz x nbsv logical array, true if the final solver attempted
%           still failed to integrate without hitting a negative, NaN, or
%           Inf value.

t = G.t;
dt = G.dt;
[nz,nx,nbsv] = size(oldbio);
nmethod = length(varargin);

count = 1;

while count <= nmethod
    
    switch varargin{count}
        
        case 'euler'
            
            [dbdt, Splitdbdt, Diag] = feval(fun, t, oldbio, G, P, B);
            newbio = oldbio + dbdt.*dt;
            
            %force biomass when necessary
            if any(B.biots(4,:)==-1) % if forcing biomasses exist
                 idx = find(B.biots(4,:) == -1); %find column indices where those forcing biomasses are
                 poolidx = B.biots(3,idx); % find corresponding critter indices
                 for i = 1:length(idx)
                     newbio(:,:,poolidx(i)) = repmat(B.biots(5,idx(i)).*G.area./(nz.*nx),...
                         [nz nx]); % divide total moleN per number of spatial boxes molN, will be properly redistributed right after
                     newbio(:,:,poolidx(i)) = newbio(:,:,poolidx(i))./(G.dz.*G.dx.*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                 end    
            end
            
            % redistribute nekton based on nekdist input
            newbio(:,:,B.isnek)= distributionnekton(newbio(:,:,B.isnek),...
                B.nekdist, G, nz, nx);
            
            failflag = false;
          
            setappdata(0, 'solvercheck', 1); % debugging
            
            
        case 'ode4'
            
            odefun = @(t,b,P) wceode(t, b, P);
            [newbio, dbdt, Splitdbdt, Diag, failflag] = ode4splitsnonneg(fun, [t t+dt], oldbio, G, P, B);
            newbio = newbio(:,:,:,2);
            
            %force biomass when necessary
            if any(B.biots(4,:)==-1) % if forcing biomasses exist
                 idx = find(B.biots(4,:) == -1); %find column indices where those forcing biomasses are
                 poolidx = B.biots(3,idx); % find corresponding critter indices
                 for i = 1:length(idx)
                     newbio(:,:,poolidx(i)) = repmat(B.biots(5,idx(i)).*G.area./(nz.*nx),...
                         [nz nx]); % divide total moleN per number of spatial boxes molN, will be properly redistributed right after
                     newbio(:,:,poolidx(i)) = newbio(:,:,poolidx(i))./(G.dz.*G.dx.*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                 end    
            end
            
            % redistribute nekton based on nekdist input
            newbio(:,:,B.isnek)= distributionnekton(newbio(:,:,B.isnek),...
                B.nekdist, G, nz, nx);

            
            setappdata(0, 'solvercheck', 2); % debugging
            
        case 'ode45'
            
            Opt = odeset('nonnegative', 1:numel(oldbio));
            [tout, newbio, dbdt, Splitdbdt] = odewrap(@ode45, fun, [t t+dt], oldbio, Opt, G, P, B);  
            failflag = false;
            
            newbio = endonly(newbio);
            
            %force biomass when necessary
            if any(B.biots(4,:)==-1) % if forcing biomasses exist
                 idx = find(B.biots(4,:) == -1); %find column indices where those forcing biomasses are
                 poolidx = B.biots(3,idx); % find corresponding critter indices
                 for i = 1:length(idx)
                     newbio(:,:,poolidx(i)) = repmat(B.biots(5,idx(i)).*G.area./(nz.*nx),...
                         [nz nx]); % divide total moleN per number of spatial boxes molN, will be properly redistributed right after
                     newbio(:,:,poolidx(i)) = newbio(:,:,poolidx(i))./(G.dz.*G.dx.*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                 end    
            end
            
            % redistribute nekton based on nekdist input
            newbio(:,:,B.isnek)= distributionnekton(newbio(:,:,B.isnek),...
                B.nekdist, G, nz, nx);

            dbdt = (newbio - oldbio)./dt;
            
            % TODO calculate proper splits by running ode4splitsnonneg at
            % each time step rather than just grabbing the first time step
            
            Splitdbdt = Splitdbdt(1);
            
            Diag = [];
            
            setappdata(0, 'solvercheck', 3); % debugging
            
        case 'implicit'
            
            [dbdt, Splitdbdt, Diag] = feval(fun, t, oldbio, G, P, B);
            
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
            lambda = dbdt./oldbio;
            lambda(oldbio == 0) = 0;
            
%             newbio = (oldbio + dbdtexplicit.*dt)./(1 - lambda.*dt);
            newbio = oldbio./(1 - lambda.*dt);
            
            %force biomass when necessary
            if any(B.biots(4,:)==-1) % if forcing biomasses exist
                 idx = find(B.biots(4,:) == -1); %find column indices where those forcing biomasses are
                 poolidx = B.biots(3,idx); % find corresponding critter indices
                 for i = 1:length(idx)
                     newbio(:,:,poolidx(i)) = repmat(B.biots(5,idx(i)).*G.area./(nz.*nx),...
                         [nz nx]); % divide total moleN per number of spatial boxes molN, will be properly redistributed right after
                     newbio(:,:,poolidx(i)) = newbio(:,:,poolidx(i))./(G.dz.*G.dx.*(G.area./sum(G.dx(1,:))));% transform into moleN.m-3
                 end    
            end
            
            % redistribute nekton based on nekdist input
            newbio(:,:,B.isnek)= distributionnekton(newbio(:,:,B.isnek),...
                B.nekdist, G, nz, nx);


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

        
    
