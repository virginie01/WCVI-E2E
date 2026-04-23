function varargout = endonly(varargin)
%ENDONLY Grab only last row of a 4-dimensional array
%
% [anew, bnew, ...] = endonly(a,b);
%
% Returns only last row of a n x m x p x q array (i.e. squeeze(a(end,:,:))).
% I wrote this function simplify grabbing the last time step result from an
% ODE solver.  Result will be a m x p x q array.
% 
% Input variables:
%
%	a,b, ...:   four-dimensional arrays.  Arrays do not have to be the
%               same size as each other.
%
% Output variables:
%
%   anew...:    three-dimensional arrays holding first row of each input
%               array, respectively

% Copyright 2009 Kelly Kearney

for ivar = 1:nargin
    temp = squeeze(varargin{ivar}(end,:,:,:));
    
    sz1 = size(temp); % nz x nx x nbsv
    sz2 = size(varargin{ivar}); % 2 x nz x nx x nbsv
    if isvector(temp) && ~isequal(sz1, sz2(2:end))
        temp = temp';
    end
    
    varargout{ivar} = temp;
    
end