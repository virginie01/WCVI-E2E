function [Ht, Wnd, VMxng, LB, Buoy, varargout] = forcingdata(Grd, it, In, varargin)
%
% This function expects all required forcing data to be stored in .mat files. 
% Each .mat file must contain a variable named after the dataset (e.g., 'QI')
% For physical forcing files:
%   - data is a numeric array where columns 1:6 are datevec (Y M D H M S)
%   - columns 7:end are the data values
%   - rows correspond to time steps
%
% It matches the current simulation time to the correct row and extracts a
% block of nsteps rows (default: 7).
%
% DIRECTORY LAYOUT:
%   Physics forcing:  fullfile(mlpath, 'data', 'physics', 'Forcing')
%     QI, AIRTMP, DEWPTT, QO, TAUY, TAUY2, WSPD10, XFIL, XFILPHY,
%     MLD, ENT, PRATE, SBC, DC, CU
%
%   Physics LB:       fullfile(mlpath, 'data', 'physics', 'LB')
%     TMP, SAL
%
%   Bio forcing:      fullfile(mlpath, 'data', 'bio', 'Forcing')
%     O2, BIOTS
%
%   Bio LB:           fullfile(mlpath, 'data', 'bio', 'LB')  (defaultdir2)
%     mixed-group boundary time series (one file per group name)
%
% OUTPUTS:
%
%   Ht:     A nested structure holding structures related to heat forcing:
%
%           QI:        incident solar radiation from Qi_input
%           AIRTMP:    air temperature from airtmp_input
%           DEWPTT:    dew point temperature from dewptT_input
%           QO:        estimate of clear sky iradiance, based on the 
%                      "Smithsonian Formula" from Seckel and Beaudry, as
%                      reported in Reed, 1977, JPO, 7, pp. 482-485. It is
%                      good for latitudes between 20S and 60N.
%
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                   x:      1 x nx array describing the boxes in the x dimension (except for Qo)
%                   data:   length(t) x length(x) array holding the data for the current time steps 
%
%
%   Wnd:    A nested structure holding structures related to wind forcing:
%           TAUY:   S-N wind stress in N.m-2
%           TAUY2:  S-N wind stress in N.m-2 from Ze
%           WSPD10: wind speed at 10m above sea level (m/s)
%           XFIL:   Remote upwelling index from Ze (m/s)
%   
%                   Each of these structures consist of the following 
%                   fields:
%
%                   t:      Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                   x:      1 x nx array describing the boxes in the x dimension 
%                   data:   length(t) x length(x) array holding the data for the current time steps 
%
%   VMxng:  A nested structure holding structures relating to vertical
%           mixing:
%
%           MLD:    Mixed layer depth coming from mld_input
%                       
%                   This structure consists of the following fields:
%                     
%                   t:      Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                   x:      1 x nx array describing the boxes in the x dimension
%                   data:   length(t) x length(x) array holding the data for the current time steps
%
%           ENT: entrainment rate from entrnmnt_input
%
%                     This structure consists of the following fields:
%                     t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                     x:    a 1 x nx array describing the boxes in the x dimension
%                     data: length(t) x length(x) array holding the data for the current time steps
%
%   LB:             structure containing boundary variables:
%
%                   TMP:   1 x 1 structure of data holding lateral boundary
%                          temperatures from LBtmp_input
%
%                          t: Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                          o: 1 x 14 array specifying the origin
%                          data: length(t) x length(o) array holding the data for the current time 
%                                steps
%
%                   SAL:   1 x 1 structure of data holding lateral boundary
%                          salinities from LBsal_input
%
%                          t: Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                          o: 1 x 14 array specifying the origin
%                          data: length(t) x length(o) array holding the data for the current time 
%                                steps
%
%
%   Buoy:           Nested structure holding structures of buyoancy fluxes
%                   or alongshore currents
%                
%                   PRATE/CU/DC/SBC
%                  
%                      t: Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%                      x: 1 x nx array describing the boxes in the x dimension (only
%                      for PRATE)
%                      data: length(t) x 1 (x length(x)) array holding the data
%
%   Bio:            Existing structure with additional fields:
%
%
%                   O2:     nz x nx numerical array containing current
%                           concentrations of O2 in mol.m-3
%
%                   BIOTS:  3 x ngroup array. Line 1 is pool id, line 2 is
%                           type id and line 3 is data for the current time
%                           step
%
%                   LB:     Nested stucture. First level = "zooplanktonic" group (i.e. mixed). 
%                           Second level:
%                           t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim 
%                                 start time)
%                           o:    1 x 14 array specifying the boundary
%                           data: length(t) x length(o) array holding the data for the current time 
%                                 steps

% -------------------------
% Time vectors
% -------------------------
tdata = Grd.time(it):In.datadt:Grd.time(it+1);

dnstart = datenum(Grd.start_date);
dncurrent = dnstart + (Grd.time(it)./86400);
dvcurrent = datevec(dncurrent);

x = {'Shelf', 'Slope'};

o_bio = {'upper open ocean','lower open ocean','rain shelf', 'rain slope',...
    'freshwater', 'vicc','sbc shelf ul','sbc shelf ll','sbc slope ul',...
            'dc shelf ul','dc shelf ll','dc slope ul','cu slope ll'};
        
o_phy = {'upper open ocean','lower open ocean','rain shelf', 'rain slope',...
    'freshwater', 'vicc ul','vicc ll','sbc shelf ul','sbc shelf ll','sbc slope ul',...
            'dc shelf ul','dc shelf ll','dc slope ul','cu slope ll'};
        
% -------------------------
% Optional Bio input
% -------------------------
if ~isempty(varargin)
    Bio = varargin{1};
end

% -------------------------
% Which datasets are needed?
% -------------------------
if ~isempty(varargin) && Bio.bioin.isnem
   datasets = {'QI', 'AIRTMP','DEWPTT','QO','TAUY','TAUY2','WSPD10','XFIL','XFILPHY','MLD',...
       'ENT','TMP','SAL','PRATE','SBC','DC','CU','O2'}; 
elseif ~isempty(varargin) && ~Bio.bioin.isnem
   datasets = {'QI', 'AIRTMP','DEWPTT','QO','TAUY','TAUY2','WSPD10','XFIL','XFILPHY','MLD',...
       'ENT','TMP','SAL','PRATE','SBC','DC','CU','O2','BIOTS'};
elseif isempty(varargin)
       datasets = {'QI', 'AIRTMP','DEWPTT','QO','TAUY','TAUY2','WSPD10','XFIL','XFILPHY',...
           'MLD','ENT','TMP','SAL','PRATE','SBC','DC','CU'}; 
end

ndatasets = length(datasets);

% -------------------------
% Grouping for outputs
% -------------------------
ht = ["QI","AIRTMP","DEWPTT","QO"];
wnd = ["TAUY","TAUY2","WSPD10","XFIL","XFILPHY"];
vmxng = ["MLD","ENT"];
lb = ["TMP","SAL"];
buoy = ["PRATE","SBC","DC","CU"];
bio = ["O2","BIOTS"];

% -------------------------
% Directory mapping
% -------------------------
mlname = mfilename('fullpath');
mlpath = fileparts(fileparts(fileparts(mlname)));

physForcingDir = fullfile(mlpath, 'data', 'physics', 'Forcing');
physLBDir      = fullfile(mlpath, 'data', 'physics', 'LB');
bioForcingDir  = fullfile(mlpath, 'data', 'bio', 'Forcing');

% Bio LB dir for mixed-group boundaries
defaultdir2    = fullfile(mlpath, 'data', 'bio', 'LB');


% -------------------------
% Main load loop
% -------------------------
for i=1:ndatasets
    
    dataset = datasets{i};

    % -------- Heat forcing --------
    if any(ht == dataset)
       Ht.(dataset).t = tdata; 
       if strcmp(dataset,"QO")==0
       Ht.(dataset).x = x;
       Ht.(dataset).data = zeros(length(tdata),length(x));
       else
       Ht.(dataset).data = zeros(length(tdata),1);       
       end        
       
       formatSpec = [dataset,'.mat'];
       file = sprintf(formatSpec);
       data = load(fullfile(physForcingDir,file));
       
       r1 = ismember(data.(dataset)(1:end-1,2:6),dvcurrent(:,2:6),'rows');
       r1 = find(r1);
       Ht.(dataset).data = data.(dataset)(r1:r1+6,7:end);
       clear data
           
    % -------- Wind forcing --------        
    elseif any(wnd == dataset)
        Wnd.(dataset).t = tdata; 
        if strcmp(dataset,"TAUY2") || strcmp(dataset,"XFIL") || strcmp(dataset,"XFILPHY")          
        Wnd.(dataset).data = zeros(length(tdata),1);
        else
        Wnd.(dataset).x = x;
        Wnd.(dataset).data = zeros(length(tdata),length(x));
        end  
        
        formatSpec = [dataset,'.mat'];
        file = sprintf(formatSpec);
        data = load(fullfile(physForcingDir,file));           
        
        r1 = ismember(data.(dataset)(1:end-1,2:6),dvcurrent(:,2:6),'rows'); 
        r1 = find(r1);
        Wnd.(dataset).data = data.(dataset)(r1:r1+6,7:end);
        clear data       
                
    % -------- Vertical mixing --------        
    elseif any(vmxng == dataset)
        VMxng.(dataset).t = tdata; 
        VMxng.(dataset).x = x;
        VMxng.(dataset).data = zeros(length(tdata),length(x));
      
        formatSpec = [dataset,'.mat'];
        file = sprintf(formatSpec);
        data = load(fullfile(physForcingDir,file));           

        r1 = ismember(data.(dataset)(1:end-1,2:6),dvcurrent(:,2:6),'rows'); 
        r1 = find(r1);
        VMxng.(dataset).data = data.(dataset)(r1:r1+6,7:end);
        clear data      
               
    % -------- Lateral boundary (physics) --------
    elseif any(lb == dataset)        
        LB.(dataset).t = tdata; 
        LB.(dataset).o = o_phy;
        LB.(dataset).data = zeros(length(tdata),length(o_phy));
        
        formatSpec = [dataset,'.mat'];
        file = sprintf(formatSpec);
        data = load(fullfile(physLBDir,file));                   

        r1 = ismember(data.(dataset)(1:end-1,2:6),dvcurrent(:,2:6),'rows');
        r1 = find(r1);
        LB.(dataset).data = data.(dataset)(r1:r1+6,7:end);
        clear data

    % -------- Buoyancy / currents --------       
    elseif any(buoy == dataset)
       Buoy.(dataset).t = tdata; 
       if strcmp(dataset,"PRATE")
       Ht.(dataset).x = x;
       Ht.(dataset).data = zeros(length(tdata),length(x));
       else
       Ht.(dataset).data = zeros(length(tdata),1);       
       end        
       
       formatSpec = [dataset,'.mat'];
       file = sprintf(formatSpec);
       data = load(fullfile(physForcingDir,file));
       
       r1 = ismember(data.(dataset)(1:end-1,2:6),dvcurrent(:,2:6),'rows');
       r1 = find(r1);
       Buoy.(dataset).data = data.(dataset)(r1:r1+6,7:end);
       clear data
        
    % -------- Bio forcing (O2, BIOTS) --------       
    elseif any(bio == dataset)
        
            formatSpec = [dataset,'.mat'];
            file = sprintf(formatSpec);
            data = load(fullfile(bioForcingDir,file));
          
            
            if strcmp(dataset, "BIOTS")
                
                t = dnstart + (cell2mat(data.(dataset).t)./86400);
                t = datevec(t);                                         
            
                r1 = ismember(t(:,2:6),dvcurrent(:,2:6),'rows');
                r1 = find(r1);
                Bio.(dataset) = [data.(dataset).col;data.(dataset).data(r1,:)];
            
            else
                t = dnstart + (data.(dataset).t./86400);
                t = datevec(t);                                         

                r1 = ismember(t(1:end-1,2:6),dvcurrent(:,2:6),'rows');
                r1 = find(r1);
                Bio.(dataset)=data.(dataset).data(:,:,r1);
            end
            
            clear data
    end

end


% -------------------------
% Mixed-group Bio LB forcing
% -------------------------
if ~isempty(varargin)
mixnames = Bio.names(Bio.ismixed,2);
nmix = length(mixnames);

 for i =1:nmix
     name = mixnames{i};
     Bio.LB.(name).t = tdata;
     Bio.LB.(name).o = o_bio;
     Bio.LB.(name).data = zeros(length(tdata),length(o_bio));
     
     if all(~strcmp(name,{'LargeZooplanktonNoDiapause','LargeZooplanktonDiapause'}))
     
     formatSpec = [name,'.mat'];
     file = sprintf(formatSpec);
     data = load(fullfile(defaultdir2,file));
     
     r1 = ismember(data.(name)(1:end-1,2:6),dvcurrent(:,2:6),'rows');
     r1 = find(r1);
     Bio.LB.(name).data = data.(name)(r1:r1+6,7:end);
     clear data
     
     end

 end
 varargout{1}=Bio;
end




end
