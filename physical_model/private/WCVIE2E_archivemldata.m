function [avg, file, ncid, vid] = WCVIE2E_archivemldata(Grd, Arch, it, data, it0)
%ARCHIVEMLDATA Write mixed_layer results to file
%
% New version should be compatible with all versions of Matlab, as long as
% snctools is installed and mexnc points to the proper backend (either
% newer internal stuff or the proper mex version).
%
% Input variables:
%
%   Grd:        1 x 1 structure holding spatial and temporal grid data for
%               physical model simulation.
%
%   Arch:       1 x 1 structure holding data related to archiving.
%
%   it:         index of current model time step.
%
%   data:       nvar x 4 cell array of data to be written to file.  Column
%               1 holds the actual data values for the current model time
%               step, which must be nz x nx x nt array.  Columns 2-4 are 
%               strings with the variables' short names, long names, and 
%               units, respectively.
%
%  it0:         index of first time step being simulated.
%
% Output variables:
%
%   avg:        nvar x 1 cell array, holding the newly-averaged values of each
%               variable.  Averages are performed over each archiving
%               period (and possibly over the geographical domain).
%
%   file:       netcdf file containing the averaged data of interest

% Copyright 2010 Kelly Kearney


% If very first step, set up average arrays

ndata= size(data,1);

if it == 1
    if Arch.space
        sz1 = cellfun(@size, data(:,1), 'uni', 0);
        avg = cellfun(@zeros, sz1, 'uni', 0);
    else
        avg=cell(ndata,1);
        avg{:}=0;
    end
else
    avg = Arch.avg;
end

% Include the current time step in the averages


for idata = 1:ndata
    if Arch.space
    avg{idata} = avg{idata} + Arch.fraction(it) .* data{idata};
    else
    avg{idata} = avg{idata} + Arch.fraction(it) .* mean(data{idata},'all');
    end
end



% Create output files

if it == it0
    
        [file, ncid, vid] = WCVIE2E_createoutputfiles(...
            Arch.outfile, ...
            Grd.z, ...
            Grd.x, ...
            Arch.dateedge(1:end-1), ...
            Arch.dateedge(2:end), ...
            Arch.middate, ...
            Arch.space,...
            data, ...
            Arch.iens);
        
else
    
    file = Arch.file;
    ncid = Arch.ncid;
    vid = Arch.vid;

end

if Arch.islast(it)
    
    isd = cellfun(@(x) isequal([Grd.nz Grd.nx], size(x)), avg);
    isde = cellfun(@(x) isequal([1 Grd.nx], size(x)), avg);
    iss = cellfun(@(x) isscalar(x), avg);
    
    start = cell(size(avg));
    [start{isd}] = deal([1 1 Arch.bin(it)]);
    [start{isde}] = deal([Arch.bin(it) 1]);
    [start{iss}] = deal(Arch.bin(it));
    
    
    for iv = 1:length(avg)
        ncwrite(file, data{iv,2}, avg{iv}, start{iv});
    end
    
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    
end



    



