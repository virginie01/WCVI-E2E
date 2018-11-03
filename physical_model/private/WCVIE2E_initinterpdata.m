function A = WCVIE2E_initinterpdata(type, data, Grd)
%INITINTERPDATA Initialize datasets for later interpolation
%
% A = initinterpdata(type, data, Grd)
%
% This function sets up forcing datasets for use for later interpolation.
% It verifies that data is in the proper format (vs time, or time and space),
% repeats climatological data if necessary, and extends data to the edges
% of the model temporal and spatial grids if necessary (using
% nearest-neighbor extrapolation).
%
% Input variables:
%
%   type:   One of the following strings:
%               'time'                 : data describes one or more variables vs time
%               'time and surface grid': data describes a single variable vs both time and
%                                        space in the x dimension only
%               'time and boundaries'  : data describes a single variable vs both time and
%                                        the origin of the data.
%               'time and full grid'   : data describes a single variable
%                                        vs both time and space in both dimensions. 
%
%   data:   Input data
%               'time'                 : n x (3+nvar) array.  Columns 1-3 hold date vectors
%                                        corresponding to each row of data, remaining
%                                        column(s) hold variable data.  If the simulation
%                                        time spans multiple years and the data only spans a
%                                        single year, the data will be treated as a
%                                        climatology and repeated for all simulation years
%
%               'time and surface grid': (n+1) x (3+nx) array. Columns 1-3 hold Year, Month, 
%                                        Day. The rest of the columns hold data for each
%                                        box in the x dimension. Line 1, columns 1-3 are left 
%                                        blank. Line 1, columns 4-end are reserved for box names 
%                                        in the x dimension. If the simulation time spans
%                                        multiple years and the data only spans a single
%                                        year, the data will be treated as a climatology and
%                                        repeated for all simulation years.
%
%
%               'time and boundaries'  :(n+1) x (3+5) array. Columns 1-3 of this array
%                                       hold a year, month and day corresponding to the dates
%                                       of the boundary conditions, and row 1 of the array
%                                       holds the origin: Column 1 = open ocean upper
%                                       layer; Column 2 = open ocean lower layer; Column 3 = 
%                                       rain; Column 4 = freshwater from run-offs; 
%                                       Column 5 = VICC. Columns 1-3 in row 1 are 
%                                       just placeholders and will be ignored. The remaining cells
%                                       hold lateral boundary data for the
%                                       given times and sources. If the simulation time spans
%                                       multiple years and the data only spans a single
%                                       year, the data will be treated as a climatology and
%                                       repeated for all simulation years.
%
%               'time and full grid'   :(nz x nx x n) x 6 cell array, with columns representing year,
%                                        month, day, the box in the x dimension, the box in the z 
%                                        dimension and the data. If the simulation time spans
%                                        multiple years and the data only spans a single
%                                        year, the data will be treated as a climatology and
%                                        repeated for all simulation years.
%
%   Grd:                                 Struct holding spatial and temporal grid data for 
%                                        WCVIE2E_physicalmodel.
%
% Output variables:
%
%   A:      1 x 1 structure with the following fields:
%
%           t:      nt x 1 vector, time (Days from simulation start
%                   time) 
%
%           z:      nz x 1, specifying the box in the z dimension ('time and full grid'
%                   only)
%
%           x:      1 x nx, specifying the box in the x dimension ('time and surface/full grid'
%                   only)
%
%           o:      origins of the boundary data
%
%           data:   nt x nvar ('time') or nt x nx ('time and surface grid') array of variables or 
%                   nz x nx x nt ('time and full grid') array of variable values 

% Copyright 2009 Kelly Kearney
% Modified by Virginie Bornarel

%--------------------------
% Setup
%--------------------------

t0 = datenum(Grd.start_date); 


switch type
    
    %--------------------------
    % Type 'time'
    %--------------------------
    
    case 'time'     
        
        % Separate data
        
        A.data = data(:,4:end);
        
        % Check for climatology
        
        dv = data(:, 1:3);
        if Grd.start_date(1) < Grd.end_date(1) && length(unique(dv(:,1))) == 1
            [dv, A.data] = repeatclimatology(dv, A.data, Grd.start_date(1), Grd.end_date(1), 1);
        end
        
        A.t = (datenum(dv) - t0);
        
        %remove the 29 of february when year is not a lap year
        [A.t,ia]=unique(A.t,'last');
        A.data=A.data(ia,:);
        
        
        % Extrapolate to edges of time grid
        
        if A.t(1) > 0
            A.t = [0; A.t];
            A.data = [A.data(1,:); A.data];
        end
        
        if A.t(end) < Grd.tmax
            A.t = [A.t; Grd.tmax];
            A.data = [A.data; A.data(end,:)];
        end
        
    %------------------------------
    % Type 'time and surface grid'
    %------------------------------
     
    case 'time and surface grid'
        
        % Separate data
        
        A.x = Grd.boxx;

        dv = cell2mat(data(2:end,1:3));
        A.t = (datenum(dv) - t0);
        A.data = cell2mat(data(2:end, 4:end));
        
        % Check for climatology
        
        if Grd.start_date(1) < Grd.end_date(1) && length(unique(dv(:,1))) == 1
            [dv, A.data] = repeatclimatology(dv, A.data, Grd.start_date(1), Grd.end_date(1), 2);
        end
        
        A.t = (datenum(dv)-t0);
        
        %remove the 29 of february when year is not a lap year
        [A.t,ia]=unique(A.t,'last');
        A.data=A.data(ia,:);
        
        % Extrapolate to edges of time grid
        
        if A.t(1) > 0
            A.t = [0; A.t];
            A.data = [A.data(1,:); A.data];
        end
        
        if A.t(end) < Grd.tmax
            A.t = [A.t; Grd.tmax];
            A.data = [A.data; A.data(end,:)];
        end

        
    %------------------------------
    % Type 'time and boundaries'
    %------------------------------
     
    case 'time and boundaries'
        
        % Separate data
        
        A.o ={'Upper open ocean';'Lower open ocean';'Rain';'Freshwater';'VICC'};

        dv = cell2mat(data(2:end,1:3));
        A.t = (datenum(dv) - t0);
        A.data = cell2mat(data(2:end, 4:end));
        
        % Check for climatology
        
        if Grd.start_date(1) < Grd.end_date(1) && length(unique(dv(:,1))) == 1
            [dv, A.data] = repeatclimatology(dv, A.data, Grd.start_date(1), Grd.end_date(1), 2);
        end
        
        A.t = (datenum(dv)-t0);
        %remove the 29 of february when year is not a lap year
        [A.t,ia]=unique(A.t,'last');
        A.data=A.data(ia,:);
        
        
        % Extrapolate to edges of time grid
        
        if A.t(1) > 0
            A.t = [0; A.t];
            A.data = [A.data(1,:); A.data];
        end
        
        if A.t(end) < Grd.tmax
            A.t = [A.t; Grd.tmax];
            A.data = [A.data; A.data(end,:)];
        end
    

    %------------------------------
    % Type 'time and full grid'
    %------------------------------
    
    case 'time and full grid'
            
        % Separate data
            
        A.z = Grd.boxz;
        A.x = Grd.boxx;
        
        dv = cell2mat(data(:,1:3));
        dvs = unique(dv,'rows');
        A.t = (datenum(dvs) - t0);
        
        A.data = zeros(length(A.z),length(A.x),length(A.t));
        
        for i = 1:size(data,1)
            
            date = (datenum(dv(i,:))-t0);
            position_date = find(A.t==date);
            
            x= data{i,4};
            position_x = find(strcmp(A.x,x));
            
            z= data{i,5};
            position_z = find(strcmp(A.z,z));
            
            A.data(position_z,position_x,position_date)=data{i,6};
        end
            
        
        % Check for climatology
        
   
        if Grd.start_date(1) < Grd.end_date(1) && length(unique(dv(:,1))) == 1
            [dv, A.data] = repeatclimatology(dvs, A.data, Grd.start_date(1), Grd.end_date(1), 3);
        end
        
            A.t = (datenum(dv)-t0);
            
        %remove the 29 of february when year is not a lap year
        [A.t,ia]=unique(A.t,'last');
        A.data=A.data(:,:,ia);
        
        
        % Extrapolate to edges of time and space grid
        
        if A.t(1) > 0
            A.t = [0; A.t];
            A.data = [A.data(:,:,1) A.data];
        end
        
        if A.t(end) < Grd.tmax
            A.t = [A.t; Grd.tmax];
            A.data = [A.data  A.data(:,:, end)];
        end

end


%--------------------------
% Repeat climatological 
% data over simulation 
% period
%--------------------------

function [dvc, datac] = repeatclimatology(dv, data, syear, eyear, flag)

% if size(data,2) == size(dv,1) && size(data,1) ~= size(dv,1)
%     flag = 1;
% elseif size(data,1) == size(dv,1)
%     flag = 2;
% else
%     error('Can''t figure out data format');
% end

yrs = syear:eyear;
nyr = length(yrs);

if flag == 1 % 'time' type
    nperyear = size(dv,1);
    dvc = repmat(dv, nyr, 1);
    dvc(:,1) = kron(yrs', ones(nperyear,1));
    datac = repmat(data, nyr, 1);
elseif flag == 2 % 'time and surface grid' or 'time and boundaries' type
    nperyear = size(dv,1);
    dvc = repmat (dv, nyr, 1);
    dvc(:,1) = kron(yrs', ones(nperyear,1));
    datac = repmat(data, nyr, 1);
else flag == 3    % 'time and full grid' type
    nperyear = size(dv,1);
    dvc = repmat(dv, nyr, 1);
    dvc(:,1) = kron(yrs', ones(nperyear,1));
    datac = repmat(data,[1 1 nyr]);

end











