function [file, ncid, vid] = diagoutput(Grd, diag, it, data, it0)
%DIAGOUTPUT writes diagnostic variables to file
%
% New version should be compatible with all versions of Matlab, as long as
% snctools is installed and mexnc points to the proper backend (either
% newer internal stuff or the proper mex version).
%
% Input variables:
%               
%   Grd:        1 x 1 structure holding spatial and temporal grid data for
%               physical model simulation 
%
%   diag:       1 x 1 structure holding current file with data up to
%               previous time step
%
%   it:         index of current model time step.
%
%   data:       nvar x 4 cell array of data to be written to file.  Column
%               1 holds the actual data values, which must be either 
%               nz x nx array, 1 x nx array or scalar. Columns 2-4 are 
%               strings with the variables'short names, long names, 
%               and units, respectively.
%
%   it0:        index of first time step being simulated.
%
% Output variables:
%
%   file:       netcdf file holding the values of each variable of interest
%               up to current time step.

% Copyright 2010 Kelly Kearney

% Create output files

if it == it0

[file, ncid, vid] = creatediagoutputfiles(...
            diag.folder,...
            Grd.z, ...
            Grd.x, ...
            Grd.time(2:end), ...
            data...
            );
        
else
    
file = diag.file;
ncid = diag.ncid;
vid = diag.vid;

end

% Write data to files
    
isd = cellfun(@(x) isequal([Grd.nz Grd.nx], size(x)), data(:,1));
isde = cellfun(@(x) isequal([1 Grd.nx], size(x)), data(:,1));
iss = cellfun(@(x) isscalar(x), data(:,1));
    
start = cell(size(data(:,1)));
[start{isd}] = deal([1 1 it]);
[start{isde}] = deal([it 1]);
[start{iss}] = deal(it);
    
    
for iv = 1:length(data(:,1))
    ncwrite(file, data{iv,2}, data{iv,1}, start{iv});
end
    


