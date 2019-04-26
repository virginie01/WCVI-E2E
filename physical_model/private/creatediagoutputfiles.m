function [file, ncid, vid] = creatediagoutputfiles(folder, z, x, t, data)

% Add dimensions
nz = size(z,1);
nx = size(x,2);
nt = length(t);


zval=-z(:,1);%normally a spatial grid would have same depths whatever the
            %location along the other dimensions. Here need to get around
            %this problem by picking either one of the 2 columns
xval=x(1,:);

Dim(1).Name = 'depth';
Dim(1).Length = nz;
Dim(1).Unlimited = false;

Dim(2).Name = 'longitude';
Dim(2).Length = nx;
Dim(2).Unlimited = false;

Dim(3).Name = 'time';
Dim(3).Length = nt;
Dim(3).Unlimited = true;

% Coordinate variable file

C.Name = '/';
C.Format = 'netcdf4';

C.Variables(1).Name = 'depth';
C.Variables(1).Dimensions = Dim(1);
C.Variables(1).Attributes = struct('Name', {'unit', 'negative'}, 'Value', {'m', 'down'});
C.Variables(1).Datatype = 'double';

C.Variables(2).Name = 'longitude';
C.Variables(2).Dimensions = Dim(2);
C.Variables(2).Attributes = struct('Name', {'unit', 'positive'}, 'Value', {'m', 'east'});
C.Variables(2).Datatype = 'double';

C.Variables(3).Name = 'time';
C.Variables(3).Dimensions = Dim(3);
C.Variables(3).Attributes = struct('Name', {'unit'}, 'Value', {'seconds since Jan-1-1992 00.00.00'});
C.Variables(3).Datatype = 'double';

% Create folder

if exist(folder, 'dir')
%     warning('Output folder already exists; existing files will not be explicitly overwritten and may cause conflicts');
else
    mkdir(folder);
end

dimfile = fullfile(folder, 'dimensions.nc');
if ~exist(dimfile, 'file')

    try
        ncwriteschema(dimfile, C);
    catch
        pause(5); % In case of weird parallel bug, pause and try again
        W = warning('off', 'all'); % b/c nc.createvariable is going to throw a bazillion warnings about variables existing
        ncwriteschema(dimfile, C);
        warning(W);
    end

    ncwrite(dimfile, 'depth', zval);
    ncwrite(dimfile, 'longitude', xval);
    ncwrite(dimfile, 'time', t);
end

% Data variables

nvar = size(data, 1);

A = struct('Name', '/', 'Format', 'netcdf4');

for iv = 1:nvar
    
    A.Variables(iv).Name = data{iv,2};
    
    if isequal(size(data{iv,1}), [nz nx])
        A.Variables(iv).Dimensions = Dim([1 2 3]);
    elseif isequal(size(data{iv,1}), [1 nx])
        A.Variables(iv).Dimensions = Dim([3 2]);
    elseif isequal(size(data{iv,1}), [1 1])
        A.Variables(iv).Dimensions = Dim(3);
    end
    
    A.Variables(iv).Attributes = struct('Name', {'long_name', 'unit'}, 'Value', data(iv,3:4));
    
    A.Variables(iv).Datatype = 'double';
    A.Variables(iv).FillValue = -9999;
    
end

file = fullfile(folder,'data_values.nc');    
if exist(file, 'file')
    warning('WCE:fileexists', 'Specified outputfile already exists; attempting to overwrite file');
    W = warning('off', 'all'); % nc.createvariable thing again
    ncwriteschema(file, A);
    warning(W);
else
    ncwriteschema(file, A);
end  
    
ncid = netcdf.open(file, 'WRITE');
vid = zeros(nvar,1);
for iv = 1:nvar
    vid(iv) = netcdf.inqVarID(ncid, data{iv,2});
end

