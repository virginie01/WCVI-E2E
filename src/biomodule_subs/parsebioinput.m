function BioIn = parsebioinput(~,~,varargin)
% PARSEBIOINPUT Parse Biological input for WCVI-E2E model
%
% This function parses biological parameters for the WCVI-E2E coupled
% physical-biogeochemical model. It supports both NEMURO-only and
% NEMURO+Ecopath model configurations. Default datasets can be overridden
% via parameter/value pairs. 
%
% OUTPUT:
% BioIn: structure containing parsed biological parameters
%
% INPUT PARAMETERS (via parameter/value pairs):
%
% --- GENERAL MODEL OPTIONS ---
% isnem : (logical) True if only NEMURO model is run [default: false]
% ecosimpp : (logical) Use Ecosim primary production [default: false]
% nomix : (logical) Turn off physical mixing [default: false]
% area : (scalar) Model surface area in m^2 [default: 2.936e10]
%
% --- INITIAL CONDITIONS ---
% bnem0 : (nz x nx x 11) Initial NEMURO state variable values in mol N.m-3
% bz0 : (nz x nx x nexz) Initial extra zooplankton values in mol N.m-3
%
% --- NEMURO CONFIGURATION ---
% NemParam : (struct) NEMURO parameter set
% odesolver : (cell) Integration methods in priority order [default: {'euler'}]
% reroute : (cell n x 5) Rerouting fluxes in NEMURO col 1:  name of flux (gpp, gra, pre, res, exx, exc, ege, 
%                               mor, dec)
%                                                   col 2:  name of original source group
%                                                   col 3:  name of original sink group
%                                                   col 4:  name of new sink group
%                                                   col 5:  fraction of flux to reroute
%
%
% --- DIAPAUSE SETTINGS ---
% diapause : (logical) Enable copepod diapause [default: false]
% dpStart  : (datenum) Downward migration start [default: July 1]
% dpEnd    : (datenum) Upward migration start [default: Jan 1]
% dpEndSpan : (int) Days for upward migration [default: 91]
% dpSpan : (int) Days to spread downward migration [default: 77]
% dpPercent : (vector) Percent of ZL transferred to the directed-swimming-ZL 
%             group *per time step* over the dpSpan period on shelf and slope[default: 10]
%
% --- FUNCTIONAL RESPONSE PARAMETERS ---
% types : (ngroup x 1 cell array of strings) Type of each group in the
%         Ecopath model
%         'n':    nekton, not affected by physical mixing 
%         'z':    zooplankton, mixed 
%         Can also be one of the Nemuro state variables that overlap with 
%         EwE ('ps', 'pl', 'zs', 'zl', 'zp', 'pon'), all of which are planktonic
%
%   x:    Vulnerability parameter used in the predator-prey functional response.
%         Accepts two formats:
%           - ngroup x 1 vector (log-transformed deviations from baseline),
%             where each entry `xi` modifies vulnerability as:
%             Xij = exp(xi + xj) + 1
%             This is analogous to the P4/P5 inputs in Aydin-Ecosim.
%             Range: -Inf to +Inf; default: 0.
%
%           - ngroup x ngroup matrix of direct (non-transformed) values
%             for each predator-prey pair Xij.
%             Range: 1 (low vulnerability) to +Inf (high vulnerability);
%             default: 2.
%
%   d:    Handling time parameter. Same accepted formats as `x`:
%           - ngroup x 1 vector of log-transformed values.
%           - ngroup x ngroup matrix of explicit handling times.
%
%   theta:Switching parameter controlling the type of functional response.
%         Same accepted formats as `x` and `d`:
%           - ngroup x 1 vector of log-transformed values, combined as:
%             THij = exp(0.5 * (thi + thj))
%
%           Functional response interpretation:
%           - THij = 1  →  Type II response (standard saturating)
%           - THij = 2  →  Type III response (sigmoid/threshold)
% grmax, thresh: (ngroup x ngroup matrix). Zooplankton maximum grazing rate (s-1) 
%                and threshold value for grazing (molN/m3).
%                To be used when the NEMURO thresholded Ivlev functional 
%                response is used instead of the planktonic foraging arena
% kgra : (vector) Temperature coefficients for extra zooplankton. Order is 
%        the same as in Ecopath group listing.
% lambda : (vector) Ivlev constants for extra zooplankton. Order is 
%          the same as in Ecopath group listing.
% m0exp : (vector/scalar) Mortality function exponent
%
% --- ECOPATH MODEL ---
% EM : (struct) Ecopath model structure


% -------------------------
% Set up input parser
% -------------------------
p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

% -------------------------
% Parse all inputs
% -------------------------

% Possibility to flag for nemuro-only
p.addParameter('isnem',        false,                             @(x) islogical(x) && isscalar(x));

% NEMURO and WCVI-E2E
p.addParameter('NemParam',     [],                                @isstruct);
p.addParameter('bnem0',        []);
p.addParameter('odesolver',    {'euler'});
p.addParameter('reroute',      cell(0,5),                         @(x) isempty(x) || (iscell(x) && size(x,2)==5));

p.addParameter('diapause',     false,                             @(x) islogical(x) && isscalar(x));
p.addParameter('dpStart',      738703,                            @(x) isnumeric(x) && isscalar(x));%July 1
p.addParameter('dpEnd',        738522,                            @(x) isnumeric(x) && isscalar(x));%Jan 1
p.addParameter('dpEndSpan',    91,                                @(x) isnumeric(x) && isscalar(x)); 
p.addParameter('dpSpan',       77,                                @(x) isnumeric(x) && isscalar(x)); 
p.addParameter('dpPercent',    10,                                @(x) isnumeric(x) && isvector(x)); 

p.addParameter('nomix',        false,                             @(x) islogical(x) && isscalar(x));
p.addParameter('types',        []);
p.addParameter('temp',         11.2967616,                        @(x) isscalar(x));
p.addParameter('kgra',         [0.084325;0.05;0.1124;0.0693;0.05],   @(x) isvector(x));
p.addParameter('lambda',       [1850;1000;2600;1350;1000],        @(x) isvector(x));
p.addParameter('ZS',           0.0000516719386931282,             @(x) isscalar(x));
p.addParameter('ZL',           3.69E-04,                          @(x) isscalar(x));
p.addParameter('m0exp',        []);
p.addParameter('EM',           [],                                @(x) validateattributes(x, {'ecopathmodel'}, {'scalar'}));
p.addParameter('x',            [],                                @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParameter('d',            [],                                @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParameter('theta',        [],                                @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParameter('grmax',        [],                                @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParameter('thresh',       [],                                @(x) validateattributes(x, {'numeric'}, {'square'}));
p.addParameter('ensdata',      [],                                @(x) validateattributes(x, {'numeric'}, {}));
p.addParameter('area',         2.936.*10^10,                      @(x) isscalar(x));       

% WCVI-E2E-only
p.addParameter('ecosimpp',     false,                             @(x) islogical(x) && isscalar(x));
p.addParameter('bz0',          []);
%p.addParameter('nekdist',     cell(0,12),                        @(x) (iscell(x) && size(x,2)==12));

% Parse
p.parse(varargin{:});
BioIn = mergestruct(p.Results, p.Unmatched);

% -------------------------
% Load datasets
% -------------------------
mlname = mfilename('fullpath');
mlpath = fileparts(fileparts(fileparts(mlname)));
paramDir = fullfile(mlpath, 'data', 'bio', 'Parameters');
icDir    = fullfile(mlpath, 'data', 'bio', 'IC');

if BioIn.isnem
    datasets = {'NemParam','bnem0','types','EM','x','d','theta'};
else
    datasets = {'NemParam','bnem0','types','EM','x','d','theta','m0exp','grmax','thresh','bz0'};
end

isdefault = ismember(datasets, p.UsingDefaults);
if ~all(isdefault) && any(isdefault)
    missingdata = sprintf('%s, ', datasets{isdefault});
    missingdata = missingdata(1:end-2);
    warning('ML:missingdata', 'Missing datasets: %s.\nYour model forcing will be a mix a your data and default data\n(which will probably lead to some weird results)', missingdata);
end

if isdefault(1)
    B = load(fullfile(paramDir, 'NemParam.mat'));  
    BioIn.NemParam = B.NemParam;
end

if isdefault(2)
    B = load(fullfile(icDir, 'bnem0.mat'));  
    BioIn.bnem0 = B.bnem0;
end

if isdefault(3)
    B = load(fullfile(paramDir, 'types.mat'));  
    BioIn.types = B.types;
end

if isdefault(4)
    B = load(fullfile(paramDir, 'EM.mat'));  
    BioIn.EM = B.EM;
end

if isdefault(5)
    B = load(fullfile(paramDir, 'x.mat'));  
    BioIn.x = B.x;
end

if isdefault(6)
    B = load(fullfile(paramDir, 'd.mat'));  
    BioIn.d = B.d;
end

if isdefault(7)
    B = load(fullfile(paramDir, 'theta.mat'));  
    BioIn.theta = B.theta;
end

if ~ BioIn.isnem
    
if isdefault(8)
    B = load(fullfile(paramDir, 'm0exp.mat'));  
    BioIn.m0exp = B.m0exp;
end

    
if isdefault(9)
    B = load(fullfile(paramDir, 'grmax.mat'));  
    BioIn.grmax = B.grmax;
end

if isdefault(10)
    B = load(fullfile(paramDir, 'thresh.mat'));  
    BioIn.thresh = B.thresh;
end

if isdefault(11)
    B = load(fullfile(icDir, 'bz0.mat'));  
    BioIn.bz0 = B.bz0;
end


end