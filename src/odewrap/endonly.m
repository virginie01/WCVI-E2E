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
%
% This file was derived from the original endonly utility developed
% by Kelly Kearney for the WCE/NEMURO framework and modified for use
% with WCVI-E2E multidimensional biological state arrays.
%
% Original framework:
% Copyright (c) 2009 Kelly Kearney
%
% Modifications by Virginie Bornarel (2022–2026) include:
%   - adaptation from 3D to 4D array handling
%   - updated indexing for WCVI-E2E state-variable structures
%   - revised dimensional documentation
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

for ivar = 1:nargin
    temp = squeeze(varargin{ivar}(end,:,:,:));
    
    sz1 = size(temp); % nz x nx x nbsv
    sz2 = size(varargin{ivar}); % 2 x nz x nx x nbsv
    if isvector(temp) && ~isequal(sz1, sz2(2:end))
        temp = temp';
    end
    
    varargout{ivar} = temp;
    
end