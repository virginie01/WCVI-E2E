function BioIn = parsewcenemin(varargin)
%PARSEWCWNEMIN Parse input for wce module
%
% Optional input variables (passed as parameter/value pairs)
%
%   isnem:      logical scalar, true if only NEMURO model is run.
%
%   NemParam:   1 x 1 structure of NEMURO parameters (see
%               nemuroParamSets.m) 
%
%   bnem0:      nz x nx x 11 array holding initial values for each NEMURO 
%               state variable (mol/m^3) represented in the model. NEMURO
%               variables which are shared with the ecopath model have to 
%               be organized in the same order as in the Ecopath model. For
%               the reamining ones, the following order holds:
%               1: PS   4: ZL   7: NH4  10: SiOH4 
%               2: PL   5: ZP   8: PON  11: Opal
%               3: ZS   6: NO3  9: DON  
%
%   bzo:        nz x nx x nexz array holding initial values for each
%               extra zooplankton group not part of NEMURO (mol/m^3).
%               Organized in the same order as in the Ecopath model.
%
%   nekdist:    nnek x 12 cell array. Each cell is a nz x nx matrix holding 
%               percentage of total biomass in each spatial box for a given 
%               nekton group (row) and month(column). nekton groups have to
%               be organized in the same order as in types/Ecopath.
%
%   odesolver:  cell array of strings, indicating which solvers to use.  If
%               the first one fails to integrate a timestep (i.e. causes
%               something to become negative, NaN, or Inf), the next one is
%               tried.  See integratebio.m for choices. [{'euler'}]
%
%   reroute:    n x 5 cell array.  This allows you to reroute fluxes from
%               the original path defined in the nemuro model.  Each row
%               desribes a change, with column as follows:
%               col 1:  name of flux (gpp, gra, pre, res, exx, exc, ege, 
%                       mor, dec)
%               col 2:  name of original source group
%               col 3:  name of original sink group
%               col 4:  name of new sink group
%               col 5:  fraction of flux to reroute
%
%   diapause:   logical scalar, true to turn on diapause for ZL group
%               [false]
%
%   
%   dpStart:    date on which diapausing copepods begin their downward
%               migration each year.  Enter as a datestr-formatted
%               month/day combo ['Sep 1']
%
%   dpEnd:      date on which diapausing copepods begin swimming back to
%               the surface.  Enter as a datestr-formatted month/day combo
%               ['Apr 1']
%    
%   dpEndSpan:  number of days over which copepods swim upward, before
%               being recombined with the non-diapausing copepod group [10]
%
%   dpSpan:     number of days over which to spread the transfer of
%               copepods from the non-diapausing group to the
%               directed-swimming group [20]
%
%  dpPercent:   Percent of passive-ZL group population to transfer to the
%               directed-swimming-ZL group *per time step* over the dpSpan
%               period. [3]
%
%  o2_input:    (nz x nx x n) x 9 array. Columns 1:6 hold year, month, day,
%               hour, minute, seconds. Columns 7 and 8 hold the name of box
%               x and z rspectively. Column 9 holds correspnding o2
%               concentration. This forcing dataset is used in the absence
%               of the 02 cycle in the model. Necessary for modeling
%               ammonififcation and denitrification
%
% H2S_input:    (nz x nx x n) x 9 array. Columns 1:6 hold year, month, day,
%               hour, minute, seconds. Columns 7 and 8 hold the name of box
%               x and z respectively. Column 9 holds correspnding H2S
%               concentration. This forcing dataset is used in the absence
%               of the H2S cycle in the model. Necessary for modeling
%               denitrification
%
% WCE only:
%
%   Ewein:      1 x 1 Ewe input structure.  Must include all ecopath fields
%               (see ecopathlite) in units of mol N/m^2/s, as well as the
%               following fields:
%                       
%               x:      Vulnerability parameter.  Can be either a ngroup x
%                       1 vector of log-transformed anomaly-from-base
%                       values (same as P4/P5 inputs in aydin-ecosim, range
%                       from -Inf to Inf w/ default of 0, xi in equation
%                       below),  or an ngroup x ngroup array of
%                       non-tranformed values for each predator-prey pair
%                       (range from 1 to Inf w/ default of 2, Xij in
%                       equation below). 
%
%                       Xij = exp(xi + xj) + 1                               
%
%               d:      Handling time parameter, same format as x. 
%
%               theta:  Switching parameter.  Same format as x but with
%                       slightly different transform  
%
%                       THij = exp(0.5 * (thi + thj))
%
%                       TH = 1 yields a type 2 functional response
%                       TH = 2 yields a type 3 functional response 
%
%   types:      ngroup x 1 cell array of strings, indicating what type of
%               critter each functional group is. 
%               'n':    nekton, not affected by physical mixing, feed over
%                       entire water column 
%               'z':    zooplankton, mixed, feed only in the layer where
%                       they are located 
%               Can also be any of the 11 nemuro state variables ('ps',
%               'pl', 'zs', 'zl', 'zp', 'no3', 'nh4', 'pon', 'don',
%               'sioh4', 'opal'), all of which are planktonic.
%
%   ecosimpp:   logical scalar, if true ecosim primary production function
%               is used, decoupling biology from nutrient constraints (for
%               debuging purposes) [false]
%
%   temp:       temperature associated with initial mass-balanced values,
%               used for grazing functional response (deg C). Annual mean 
%               temperature corresponding to starting year. [0] 
%
%   kgra:       nexz x 1 array, temperature coefficient associated with
%               non-nemuro zooplankton groups (i.e. those with type 'z').
%               Default is 0.0693 deg C^-1 for all, i.e. a Q10 of 2. Order
%               is the same as in Ecopath group listing.
%
%   m0exp:      Exponent for mortality function, of form M0 = aB^(m0exp). A
%               value of 1 leads to linear mortality and 2 to quadratic
%               mortality. Can also be a nlive x 1 array of values to allow
%               different functions for each critter. Ordered in the same
%               order as in Ecopath.[2] 
%
%   area:       total area of Ecopath model (m^2).
%
% Copyright 2014 Kelly Kearney

% Set up

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

% Flag for nemuro-only (as opposed to water column ecosystem)

p.addParamValue('isnem',        true,                      @(x) islogical(x) && isscalar(x));

% Both

p.addParamValue('NemParam',     [],                        @isstruct);
p.addParamValue('bnem0',        []);
p.addParamValue('bz0',          []);
p.addParamValue('nekdist',      cell(0,12),                @(x) (iscell(x) && size(x,2)==12));
p.addParamValue('odesolver',    {'euler'});
p.addParamValue('reroute',      cell(0,5),                 @(x) isempty(x) || (iscell(x) && size(x,2)==5));


p.addParamValue('diapause',     false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('dpStart',      'Sep 1',                   @(x) ischar(x));
p.addParamValue('dpEnd',        'Apr 1',                   @(x) ischar(x));
p.addParamValue('dpEndSpan',    10,                        @(x) isnumeric(x) && isscalar(x)); 
p.addParamValue('dpSpan',       20,                        @(x) isnumeric(x) && isscalar(x)); 
p.addParamValue('dpPercent',    3,                         @(x) isnumeric(x) && isscalar(x)); 

% Wce-only

p.addParamValue('nomix',        false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('ecosimpp',     false,                     @(x) islogical(x) && isscalar(x));
% p.addParamValue('Ewein',        [],                        @isstruct);
p.addParamValue('types',        []);
p.addParamValue('temp',         0,                         @(x) isscalar(x));
p.addParamValue('kgra',         0.0693,                    @(x) isvector(x));
p.addParamValue('m0exp',        2,                         @(x) isnumeric(x) && (isscalar(x) || isvector(x)));
% p.addParamValue('predatdepth',  true,                      @(x) islogical(x) && isscalar(x));
p.addParamValue('EM',           [],                        @(x) validateattributes(x, {'ecopathmodel'}, {'scalar'}));
p.addParamValue('x',            [],                        @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParamValue('d',            [],                        @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParamValue('theta',        [],                        @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParamValue('ensdata',      [],                        @(x) validateattributes(x, {'numeric'}, {}));
p.addParamValue('area',         3.1.*10^10,                @(x) isscalar(x));       

% requried biogeochemical datasets 

p.addParamValue('o2_input', [], @(x) size(x,2) == 9);
p.addParamValue('H2S_input', [], @(x) size(x,2) == 9);


% Parse

p.parse(varargin{:});
BioIn = mergestruct(p.Results, p.Unmatched);

% Add default datasets if necessary

mlname = mfilename('fullpath');
mlpath = fileparts(fileparts(mlname));
defaultdir = fullfile(mlpath, 'defaultdata');

datasets = {'o2_input','H2S_input'};
isdefault = ismember(datasets, p.UsingDefaults);
if ~all(isdefault) && any(isdefault)
    missingdata = sprintf('%s, ', datasets{isdefault});
    missingdata = missingdata(1:end-2);
    warning('ML:missingdata', 'Missing datasets: %s.\nYour model forcing will be a mix a your data and default data\n(which will probably lead to some weird results)', missingdata);
end

if isdefault(1)
    B = load(fullfile(defaultdir, 'o2_input.mat'));  
    BioIn.o2_input = B.o2input;
end

if isdefault(2)
    B = load(fullfile(defaultdir, 'H2S_input.mat'));   
    BioIn.H2S_input = B.H2Sinput;
end

% Check that required inputs are here

if BioIn.isnem
    nodefault = {'NemParam','bnem0'};
else
%     nodefault = {'Ewein', 'NemParam', 'types', 'bnem0'};
    nodefault = {'EM', 'NemParam', 'types', 'bnem0','bz0','nekdist', 'x', 'd', 'theta'};
end

tf = ismember(p.UsingDefaults, nodefault);
missing = sprintf('%s, ', p.UsingDefaults{tf});
if any(tf)
    error('Missing input for wce/nemuro module: %s', missing(1:end-2));
end
