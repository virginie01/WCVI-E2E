function In = WCVIE2E_parseinput(varargin)   
%PARSEINPUT Parse and check all input variables
%
% In = parseinput(varargin) 
%
% This function checks all input variables passed to the WCVIE2E_physicalmodel.
% It supplies default values where applicable, and runs some simple
% validation checks on the values of each variable. See WCVIE2E_physicalmodel
% for input variable list.
% 
% Input variables:
%
%	see WCVIE2E_physicalmodel
%
% Output variables:
%
% 	In:	structure holding all input variables, with default values substituted
% 		where applicable   

% Copyright 2012 Kelly Kearney 		


p = inputParser;

isnumscal = @(x) isnumeric(x) && isscalar(x);

% No partial matches; it can cause a seg fault w/ some of my options

if isprop(p, 'PartialMatching') % only R2013+
    p.PartialMatching = false;
end

% The only required input: the output file name

p.addRequired('outputfile', @ischar);

% Model parameters and their defaults

p.addParamValue('dz', [25 42.5; 85 347.5; 10 10], @(x) isnumeric(x) && ismatrix(x));
p.addParamValue('dx', [20000 10000; 20000 10000; 20000 10000] , @(x) isnumeric(x) && ismatrix(x));
p.addParamValue('dt', 1, isnumscal);
p.addParamValue('syear', 1992, @(x) isscalar(x) || isequal(size(x), [1 6]));
p.addParamValue('eyear', 2017, @(x) isscalar(x) || isequal(size(x), [1 6]));
p.addParamValue('krad1', 0.15, isnumscal);
p.addParamValue('prad1', 0.45, isnumscal);
p.addParamValue('krad2', 1.67, isnumscal);
p.addParamValue('alb', 0.079, isnumscal);
p.addParamValue('Lat', 45, isnumscal);
p.addParamValue('whgt', 10, isnumscal);
p.addParamValue('tarch', 1*30, @(x) (isnumeric(x) && isscalar(x)) || (isnumeric(x) && isvector(x)));
p.addParamValue('srelaxtime', 1*30, isnumscal);
p.addParamValue('trelaxtime', 1*30, isnumscal);
p.addParamValue('brelaxtime', 1*30, @isvector);
p.addParamValue('openbottom', false, @(x) isscalar(x) && islogical(x));
p.addParamValue('beginarchive', NaN, @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParamValue('endarchive',   NaN, @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParamValue('outputextension', [], @(x) isempty(x) || (iscell(x) && all(cellfun(@ischar,x)))); 
p.addParamValue('verbose', true, @(x) isscalar(x) && islogical(x));
p.addParamValue('hotstartdn', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
p.addParamValue('hotstartload', '',@(x) ischar(x));
p.addParamValue('hotstartsave', 'mlhotstart', @(x) ischar(x));
p.addParamValue('iens', 1, isnumscal);
p.addParamValue('nens', 1, isnumscal);
% p.addParamValue('newfile', true, @(x) isscalar(x) && islogical(x));

% Optional model parameters (default to empty)

p.addParamValue('biofun',  [], @(x) isempty(x) || isa(x, 'function_handle'));
p.addParamValue('srelax',  [], @(x) isnumeric(x) && ismatrix(x) && size(x,2) == 6);
p.addParamValue('trelax',  [], @(x) isnumeric(x) && ismatrix(x) && size(x,2) == 6);

% Placeholders for forcing datasets

p.addParamValue('Qi_input', [], @(x) size(x,2) == 5);
p.addParamValue('airtmp_input', [], @(x) size(x,2) == 5);
p.addParamValue('dewptT_input', [], @(x) size(x,2) == 5);
p.addParamValue('uWndSpd_input', [], @(x) size(x,2) == 5);
p.addParamValue('vWndSpd_input', [], @(x) size(x,2) == 5);
p.addParamValue('wndl_input', [], @(x) size (x,2) == 4);
p.addParamValue('wndr_input', [], @(x) size(x) == 4);
p.addParamValue('mld_input', [], @(x) size(x,2) == 5);
p.addParamValue('entrnmnt_input', [], @(x) size(x,2) == 6);

p.addParamValue('t_input', [], @isnumeric);
p.addParamValue('s_input', [], @isnumeric);

p.addParamValue('LBtmp_input', [], @(x) size(x,2) == 5);
p.addParamValue('LBsal_input', [], @(x) size(x,2) == 5);

% Allow structure input

p.StructExpand = true;

% Allow unmatched inputs (for biological module inputs, which could be
% anything)

p.KeepUnmatched = true;

% Parse input

p.parse(varargin{:});
In = mergestruct(p.Results, p.Unmatched);

% Add default datasets if necessary

mlname = mfilename('fullpath');
mlpath = fileparts(fileparts(mlname));
defaultdir = fullfile(mlpath, 'defaultdata');

datasets = {'Qi_input','airtmp_input','dewptT_input','uWndSpd_input',...
    'vWndSpd_input', 'wndl_input', 'wndr_input','mld_input',...
    'entrnmnt_input', 't_input', 's_input','LBtmp_input','LBsal_input'};
isdefault = ismember(datasets, p.UsingDefaults);
if ~all(isdefault) && any(isdefault)
    missingdata = sprintf('%s, ', datasets{isdefault});
    missingdata = missingdata(1:end-2);
    warning('ML:missingdata', 'Missing datasets: %s.\nYour model forcing will be a mix a your data and default data\n(which will probably lead to some weird results)', missingdata);
end

if isdefault(1)
    B = load(fullfile(defaultdir, 'Qi_input.mat'));  % loads n x 5 matrix of mld data, see above
    In.Qi_input = B.Qiinput;
end

if isdefault(2)
    B = load(fullfile(defaultdir, 'airtmp_input.mat'));  % loads n x 5 matrix of mld data, see above
    In.airtmp_input = B.airtmpinput;
end

if isdefault(3)
    B = load(fullfile(defaultdir, 'dewptT_input.mat'));  % loads n x 5 matrix of mld data, see above
    In.dewptT_input = B.dewptTinput;
end

if isdefault(4)
    B = load(fullfile(defaultdir, 'uWndSpd_input.mat'));  % loads n x 5 matrix of mld data, see above
    In.uWndSpd_input = B.uWndSpdinput;
end

if isdefault(5)
    B = load(fullfile(defaultdir, 'vWndSpd_input.mat'));  % loads n x 5 matrix of mld data, see above
    In.vWndSpd_input = B.vWndSpdinput;
end

if isdefault(6)
    B = load(fullfile(defaultdir, 'wndl_input.mat')); % loads n x 6 matrix of entrainment rate data, see above
    In.wndl_input = B.wndlinput;
end

if isdefault(7)
    B = load(fullfile(defaultdir, 'wndr_input.mat'));   % loads n x 4 matrix of remote wind forcing, see above
    In.wndr_input = B.wndrinput;
end

if isdefault(8)
    B = load(fullfile(defaultdir, 'mld_input.mat'));   % loads 3 x 2 matrix of initial temperature data, see above
    In.mld_input = B.mldinput;
end

if isdefault(9)
    B = load(fullfile(defaultdir, 'entrnmnt_input.mat'));   % loads 3 x 2 matrix of initial salinity data, see above
    In.entrnmnt_input = B.entrnmntinput;
end

if isdefault(10)
    B = load(fullfile(defaultdir, 't_input.mat'));   % loads 3 x 2 matrix of initial salinity data, see above
    In.t_input = B.tinput;
end

if isdefault(11)
    B = load(fullfile(defaultdir, 's_input.mat'));   % loads 3 x 2 matrix of initial salinity data, see above
    In.s_input = B.sinput;
end

if isdefault(12)
    B = load(fullfile(defaultdir, 'LBtmp_input.mat'));   % loads 3 x 2 matrix of initial salinity data, see above
    In.LBtmp_input = B.LBtmpinput;
end

if isdefault(13)
    B = load(fullfile(defaultdir, 'LBsal_input.mat'));   % loads 3 x 2 matrix of initial salinity data, see above
    In.LBsal_input = B.LBsalinput;
end

% Allow .nc extension to be left off of output file name (this is mainly for 
% back-compatibility with old single-file output... the .nc extension is 
% removed when I set up the archiving folder in initialize.m)

[blah,blah,ext] = fileparts(In.outputfile);
if isempty(ext) || ~strcmp(ext, 'nc')
    In.outputfile = [In.outputfile '.nc'];
end

% Indicators of optional parameters

In.hasbio     = ~isempty(In.biofun);
In.hassrelax  = ~isempty(In.srelax);
In.hastrelax  = ~isempty(In.trelax);

% Check for agreement between archiving variables

In.tarch        = In.tarch(:);
In.beginarchive = In.beginarchive(:);
In.endarchive   = In.endarchive(:);

len = cellfun(@length, {In.tarch, In.beginarchive, In.endarchive});
maxlen = max(len);
if ~all(len == maxlen | len == 1)
    error('Archiving variables must agree in number');
end
if maxlen > 1
    if len(1) == 1
        In.tarch = ones(maxlen,1)*In.tarch;
    end
    if len(2) == 1
        In.beginarchive = ones(maxlen,1)*In.beginarchive;
    end
    if len(3) == 1
        In.endarchive = ones(maxlen,1)*In.endarchive;
    end
end

if isempty(In.outputextension)
    In.outputextension = cellstr(num2str((1:length(In.tarch))'));
else
    if length(In.outputextension) ~= length(In.tarch)
        error('Incorrect number of output extensions');
    end
end

% Hot start parameters (save or use)

if ~isempty(In.hotstartload) && ~exist(In.hotstartload, 'file')
    error('Hot start file not found');
end
    
    
    