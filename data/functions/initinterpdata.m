function A = initinterpdata(type, data, Grd)
%WCVIE2E_INITINTERPDATA Prepare forcing datasets for interpolation
%
% A = WCVIE2E_initinterpdata(type, data, Grd)
%
% This function formats external forcing datasets so they can be
% interpolated consistently during the WCVI-E2E simulation.
%
% It performs the following tasks:
%   - Parses calendar dates and converts them to model time (seconds)
%   - Detects and repeats climatological forcing over the entire
%   simulation period
%   - Extends forcing data to the temporal bounds of the model grid
%     using nearest-neighbor extrapolation
%
% -------------------------------------------------------------------------
% FORCING TYPE DEFINITIONS
% -------------------------------------------------------------------------
%
% TYPE = 'time'
%   Data varies with time only
%   DATA format: n x (6 + nvar) cell array
%     Columns 1–6 : [year month day hour minute second]
%     Columns 7+  : variable values
%
% TYPE = 'time and surface grid'
%   Data varies with time and x-dimension only
%   DATA format: (n+1) x (6 + nx) cell array
%     Row 1, cols 7+ : x-box names
%     Rows 2+, cols 1–6 : date vector
%     Rows 2+, cols 7+ : data values
%
% TYPE = 'time and boundaries'
%   Data varies with time and boundary source
%   DATA format: (n+1) x 12 cell array
%     Columns 1–6 : date vector
%     Columns 7–20: boundary forcing values
%
% TYPE = 'time and full grid'
%   Data varies with time, x, and z
%   DATA format: (nz * nx * n) x 9 cell array
%     Columns 1–6 : date vector
%     Column 7    : x-box name
%     Column 8    : z-box name
%     Column 9    : data value
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% type : string
% data : cell array
% Grd  : grid structure from physicalmodel
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% A : struct with fields
%     t     : time (seconds since simulation start)
%     data  : formatted forcing data
%     x, z  : spatial labels (if applicable)
%     o     : boundary labels (if applicable)
%
% -------------------------------------------------------------------------
% NOTE
% -------------------------------------------------------------------------
% Each helper function below is a direct extraction of the corresponding
% switch-case block from the original implementation. No numerical logic,
% ordering, or behavior has been changed.
%
% Original implementation: Kelly Kearney (2009)
% WCVI-E2E adaptation: Virginie Bornarel



% -------------------------------------------------------------------------
% Validate forcing type
% -------------------------------------------------------------------------
validTypes = {
    'time'
    'time and surface grid'
    'time and boundaries'
    'time and full grid'
};

assert(any(strcmp(type, validTypes)), ...
    'WCVIE2E_initinterpdata:InvalidType', ...
    'Unknown forcing type: %s', type);

t0 = datenum(Grd.start_date);

% -------------------------------------------------------------------------
% Dispatch by forcing type
% -------------------------------------------------------------------------
switch type

    case 'time'
        A = parseTimeOnly(data, Grd, t0);

    case 'time and surface grid'
        A = parseTimeSurfaceGrid(data, Grd, t0);

    case 'time and boundaries'
        A = parseTimeBoundaries(data, Grd, t0);

    case 'time and full grid'
        A = parseTimeFullGrid(data, Grd, t0);
end

end

function A = parseTimeOnly(data, Grd, t0)

dv   = cell2mat(data(2:end,1:6));
vals = cell2mat(data(2:end,7:end));

[dv, vals] = handleClimatology(dv, vals, Grd, 1);

A.t    = (datenum(dv) - t0) * 86400;
A.data = vals;

[A.t, A.data] = extendTimeBounds(A.t, A.data, Grd.tmax);

end

function A = parseTimeSurfaceGrid(data, Grd, t0)

dv   = cell2mat(data(2:end,1:6));
vals = cell2mat(data(2:end,7:end));

[dv, vals] = handleClimatology(dv, vals, Grd, 2);

A.t    = (datenum(dv) - t0) * 86400;
A.x    = Grd.boxx;
A.data = vals;

[A.t, A.data] = extendTimeBounds(A.t, A.data, Grd.tmax);

end

function A = parseTimeBoundaries(data, Grd, t0)

A.o = {
    'Upper open ocean'
    'Lower open ocean'
    'Rain shelf'
    'Rain slope'
    'Freshwater'
    'VICC shelf UL'
    'VICC shelf LL'
    'SBC shelf UL'
    'SBC shelf LL'
    'SBC slope UL'
    'DC shelf UL'
    'DC shelf LL'
    'DC slope UL'
    'CU slope LL'
};

dv   = cell2mat(data(2:end,1:6));
vals = cell2mat(data(2:end,7:end));

[dv, vals] = handleClimatology(dv, vals, Grd, 2);

A.t    = (datenum(dv) - t0) * 86400;
A.data = vals;

[A.t, A.data] = extendTimeBounds(A.t, A.data, Grd.tmax);

end

function A = parseTimeFullGrid(data, Grd, t0)

A.x = Grd.boxx;
A.z = Grd.boxz;

dvAll = cell2mat(data(:,1:6));
dv    = unique(dvAll, 'rows');

A.t = (datenum(dv) - t0) .* 86400;
A.data = zeros(length(A.z), length(A.x), length(A.t));

for i = 1:size(data,1)
    tidx = find(A.t == (datenum(dvAll(i,:)) - t0) * 86400);
    xidx = find(strcmp(A.x, data{i,7}));
    zidx = find(strcmp(A.z, data{i,8}));
    A.data(zidx, xidx, tidx) = data{i,9};
end

[dv, A.data] = handleClimatology(dv, A.data, Grd, 3);

[A.t, A.data] = extendTimeBounds(A.t, A.data, Grd.tmax);

end

function [dvOut, dataOut] = handleClimatology(dv, data, Grd, mode)

if Grd.start_date(1) < Grd.end_date(1) && numel(unique(dv(:,1))) == 1
    yrs = Grd.start_date(1):Grd.end_date(1);
    nyr = numel(yrs);
    nper = size(dv,1);

    dvOut = repmat(dv, nyr, 1);
    dvOut(:,1) = kron(yrs', ones(nper,1));

    switch mode
        case {1,2}
            dataOut = repmat(data, nyr, 1);
        case 3
            dataOut = repmat(data, [1 1 nyr]);
    end
else
    dvOut   = dv;
    dataOut = data;
end

end

function [t, data] = extendTimeBounds(t, data, tmax)

if t(1) > 0
    t = [0; t];
    if ndims(data) == 2
        data = [data(1,:); data];
    else
        data = [data(:,:,1) data];
    end
end

if t(end) < tmax
    t = [t; tmax];
    if ndims(data) == 2
        data = [data; data(end,:)];
    else
        data = cat(3, data, data(:,:,end));
    end
end

end