function varargout = biomodel(action, varargin)
%BIOMODEL biological module
%
% This module simulates a mixed planktonic-nektonic ecosystem, based on
% a combination of models derived mainly from NEMURO and Kerim Aydin's
% version of Ecosim, with a little bit of COBALT thrown in for flavor.
%
% A note on units: The biomass of all critters is saved to file in mol
% N[Si][Fe]/m^3.  For nektonic critters, all biomass is placed in the
% surface cell, and actually represents the total over the entire water
% column; multiply by the thickness of the surface layer to get the true
% biomass, in mol N/m^2.
%
% See biomodule.m for function syntax descriptions.  For fields that must be 
% present in the In structure (passed to mixed_layer as parameter/value pairs), 
% see the help for parsewcenemin.m. true false true false false

% Copyright 2008-2015 Kelly Kearney

%nin(1) = nargin(@init) + 1;
%nin(2) = nargin(@sourcesink) + 1;
%nin(3) = nargin(@vertmove) + 1;

nout(1) = nargout(@init);
nout(2) = nargout(@sourcesink);
nout(3) = nargout(@vertmove);

switch action
    case 'init'
        
        out = cell(1, nout(1));
%       narginchk(nin(1), nin(1));      
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
%       narginchk(nin(2), nin(2));       
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
%       narginchk(nin(3), nin(3));      
        [out{:}] = vertmove(varargin{:});
        
    otherwise
        
        error('Invalid action for biological module');
end

[varargout{1:nargout}] = out{:};
        
%**************************************************************************

function [BioIn, bio, ismixed, Biovars, names, diagnames] = init(In, Grd)

%------------------------------
% Parse and check input
%------------------------------

BioIn = parsebioinput(Grd, In);

%------------------------------
% Variable names
%------------------------------

[names, nbsv, Biovars] = setstatevars(BioIn); 

%----------------------------
% Indicators for mixing and 
% bottom-forcing
%----------------------------

% Plankton mixed, nekton not
% If no-mix (debugging), nothing is mixed.
if BioIn.isnem && ~BioIn.nomix
    ismixed = true(nbsv,1);
elseif BioIn.isnem && BioIn.nomix
    ismixed = false(nbsv,1);
elseif ~ BioIn.isnem && BioIn.nomix
    ismixed = false(nbsv,1);
elseif ~BioIn.isnem && ~BioIn.nomix
    ismixed = ~Biovars.isnek;
end

if Biovars.diapause
    ismixed([Biovars.idx.zl1 Biovars.idx.zl2]) = true;
end

%---------------------------
% diagnostic variables
%---------------------------

[diagnames, Biovars] = setdiagnosticsvars(BioIn, names, ismixed, Biovars);

%----------------------------
% Initial biomass for
% state variables
%----------------------------

bio = zeros(Grd.nz, Grd.nx, nbsv);

% Nemuro-derived variables 

bio(:,:,Biovars.nemidx(1:11)) = BioIn.bnem0(:,:,:); % mol/m^3

% Extra zooplankton groups 
if ~BioIn.isnem
bio(:,:,Biovars.extrazooidx) = BioIn.bz0(:,:,:);  % mol N/m^3
end
% nekton groups

EM = BioIn.EM;
if isempty(BioIn.ensdata)
    Ep = EM.ecopath;
else
    [~, Ep] = EM.ecopath('ensemble', BioIn.ensdata);
end

% change units
Ep.b = Ep.b.*0.001885; % from t.km-2 into molesN.m-2
Ep.q0 = (Ep.q0.*0.001885)./31536000; %from t.km-2.yr-1 into molesN.m-2.s-1
Ep.pb = Ep.pb./31536000;
Ep.otherMortRate = Ep.otherMortRate./31536000;

if ~ BioIn.isnem    
bnek = Ep.b(Biovars.isnek).*BioIn.area; % mol N.m-2 * m^2= molN
%Nekton biomass stored in demersal layer where volume doesn't change over
%time. Although actually per area, here stored as per volume for
%consistency
bnek = bnek./(In.dz(3,1).*In.dx(3,1).*(BioIn.area./(In.dx(1,1)+In.dx(1,2))));%mol N.m-3-
bio(3,1,Biovars.isnek) = reshape(bnek,1,1,length(bnek));
end
% PL silica is proportional to PL N

pln = bio(:,:,Biovars.idx.pl);
plsi = pln .* BioIn.NemParam.RSiN;
bio(:,:,Biovars.idx.plsi) = plsi;

% If diapause, split ZL biomass

if Biovars.diapause
    bio(:,:, Biovars.idx.zl1) = bio(:,:, Biovars.idx.zl);
    bio(:,:, Biovars.idx.zl) = 0;
end

%----------------------------
% Variables needed for ODE
%----------------------------

% Params shared with nemurokak

[Biovars, Np] = setbioparams(BioIn, Biovars.nemidx, nbsv, Grd, Biovars);

if BioIn.isnem
    % number of groups and their idx: groups in NEMURO that overlap with EwE model
    [isEwEvar,EwEidx] = ismember(lower(names(:,1)),BioIn.types);
    EwEidx = EwEidx(isEwEvar);
    % number of groups and their idx: living organisms in NEMURO that overlap with EwE model
    [isEwEliving,EwElivingidx] = ismember({'ps';'pl';'zs';'zl';'zp'},BioIn.types);
    EwElivingidx = EwElivingidx(isEwEliving);
    isEwEliving = [isEwEliving;false(nbsv-5,1)];
end


% Mortality exponent
if BioIn.isnem
    Biovars.m0exp              = zeros(nbsv,1);
    Biovars.m0exp(isEwEliving) = 2;             % no unit
else
    Biovars.m0exp             = zeros(nbsv,1);
    Biovars.m0exp(1:EM.nlive) = BioIn.m0exp;    % no unit
end

% Nekton: the parameters for nekton come from the Ecopath mass balance
if BioIn.isnem
   Biovars.b0           = zeros(nbsv,1);
   Biovars.b0(isEwEvar) = Ep.b(EwEidx); % molN.m^-2.     
else
   Biovars.b0              = zeros(nbsv,1);
   Biovars.b0(1:EM.ngroup) = Ep.b; % molN.m^-2. 
end

if BioIn.isnem
    Biovars.q0 = zeros(nbsv);
    Biovars.q0(isEwEliving,isEwEliving) = Ep.q0(EwElivingidx,EwElivingidx);
    Biovars.ge = zeros(nbsv,1);
    Biovars.ge(isEwEliving) = Ep.ge(EwElivingidx);
    Biovars.gs = zeros(nbsv,1);
    Biovars.gs(isEwEliving) = EM.groupdata.gs(EwElivingidx);
else
   Biovars.q0 = zeros(nbsv);
   Biovars.q0(1:EM.ngroup,1:EM.ngroup) = Ep.q0(1:EM.ngroup,1:EM.ngroup);
   Biovars.import = zeros(nbsv,1);
   Biovars.import(1:EM.ngroup) = (EM.groupdata.qb.*EM.groupdata.import)./31536000; % s-1
   Biovars.ge = zeros(nbsv,1);
   Biovars.ge(1:EM.nlive) = Ep.ge(1:EM.nlive);
   Biovars.gs = zeros(nbsv,1);
   Biovars.gs(1:EM.nlive) = EM.groupdata.gs(1:EM.nlive);
end

% For x, d, and theta, I accept data in one of two ways:
% 1) vector of P-value log-tranformed-anomaly-from-base values, identical
% to those used in the functional response file input for aydin-ecosim
% 2) matrix of values for each group pair.  These values are NOT anomalies
% but actual values to be used in the functional response equations.

dcmask = table2array(EM.dc) == 0;

if isvector(BioIn.x)
    [xj,xi] = meshgrid(BioIn.x);
    x = 1 + exp(xi + xj);
    x(dcmask) = 0;
else
    x = BioIn.x;
end

if isvector(BioIn.d)
    [dj,di] = meshgrid(BioIn.d);
    d = 1 + exp(di + dj);
    d(dcmask) = 0;
else
    d = BioIn.d;
end

if isvector(BioIn.theta)
    [thj,thi] = meshgrid(BioIn.theta);
    theta = exp(0.05*(thi + thj)); % TODO double-check theta, since different in spreadsheet and ppt
    theta(dcmask) = 0;
else
    theta = BioIn.theta;
end

if BioIn.isnem
    Biovars.x = zeros(nbsv);
    Biovars.x(isEwEvar,isEwEvar)=x(EwEidx,EwEidx);
    Biovars.d = zeros(nbsv);
    Biovars.d(isEwEvar,isEwEvar)=d(EwEidx,EwEidx);
    Biovars.theta = zeros(nbsv);
    Biovars.theta(isEwEvar,isEwEvar)=theta(EwEidx,EwEidx);
    
else
    n = nbsv - EM.ngroup;
    Biovars.x = padarray(x, [n n], 'post');
    Biovars.d = padarray(d, [n n], 'post');
    Biovars.theta = padarray(theta, [n n], 'post');
end

% Add ngroup x ngroup array grmax and thresh in case I don't use the
% planktonic foraging arena

if ~BioIn.isnem
grmax = BioIn.grmax;
thresh = BioIn.thresh;
n = nbsv - EM.ngroup;
Biovars.grmax = padarray(grmax, [n n], 'post');
Biovars.thresh = padarray(thresh, [n n], 'post');
end

% Zooplankton: grazing same as predation, but in volumetric terms

Biovars.b0v = Biovars.b0./560;% mol N m^-3 based on domain area 2.936.*10^10 m2 and bottom depths 120m+1000m. But I don't think this figure matter since the curve tries to reproduce the 
Biovars.q0v = Biovars.q0./560; % mol N m^-3 s^-1

% adjust b0v and q0v for LTL (i.e. concentrations and fluxes taken from a
% one box in conversions.xlsx- ULsh June)

Biovars.b0v([1:10,55],1)=[0.00108663367294276;0.00217203605898552;0.0000516719386931282;
    0.000368604351936326;0.000191354202990411;2.72316659291623E-09;0.0000133490872197612;
    9.06785338173271E-06;5.1403304025472E-08;0.00023311108595052;0.00019060631829999];

Biovars.q0v(1,3)=2.54346236617204E-10;Biovars.q0v(1,4)=6.86053945887653E-10;
Biovars.q0v(2,4)=2.72102101053744E-09;Biovars.q0v(3,4)=6.89678646626627E-11;
Biovars.q0v(2,5)=1.9744366271254E-10;Biovars.q0v(3,5)=1.07893398643476E-11;
Biovars.q0v(4,5)=3.59258709732891E-10;Biovars.q0v(1,6)=1.10E-13;
Biovars.q0v(2,6)=5.93E-14;Biovars.q0v(3,6)=4.56E-15;Biovars.q0v(55,6)=2.41E-15;
Biovars.q0v(1,7)=2.34E-11;Biovars.q0v(2,7)=3.45E-11;Biovars.q0v(3,7)=4.83E-12;
Biovars.q0v(4,7)=3.17E-11;Biovars.q0v(5,7)=7.17E-12;Biovars.q0v(10,7)=6.66E-12;
Biovars.q0v(55,7)=9.52E-12;Biovars.q0v(1,8)=1.30E-11;Biovars.q0v(2,8)=1.44E-11;
Biovars.q0v(3,8)=2.09E-12;Biovars.q0v(4,8)=1.39E-11;Biovars.q0v(5,8)=5.01E-12;
Biovars.q0v(6,8)=0.00E+00;Biovars.q0v(7,8)=0.00E+00;Biovars.q0v(8,8)=0.00E+00;
Biovars.q0v(9,8)=0.00E+00;Biovars.q0v(10,8)=1.85E-11;Biovars.q0v(55,8)=2.26E-11;
Biovars.q0v(1,9)=3.55E-13;Biovars.q0v(2,9)=4.30E-13;Biovars.q0v(3,9)=1.86E-14;
Biovars.q0v(4,9)=1.62E-13;Biovars.q0v(55,9)=8.18E-14;Biovars.q0v(1,10)=9.10E-10;
Biovars.q0v(2,10)=6.98E-10;Biovars.q0v(3,10)=4.86E-11;
Biovars.q0v(4,10)=3.40E-10;Biovars.q0v(10,10)=2.01E-10;

% Production (for debugging with ecosim production only, adjust p0v as above when running ecosimpp)
if BioIn.isnem
    Biovars.p0 = zeros(nbsv,1);
    Biovars.p0(isEwEvar)= Ep.pb(EwEidx).*Ep.b(EwEidx);% mol N m^-2 s^-1
    Biovars.p0v = Biovars.p0./560;% mol N m^-3 s^-1    
else
    Biovars.p0 = zeros(nbsv,1);
    Biovars.p0(1:EM.ngroup) = Ep.pb .* Ep.b; % mol N m^-2 s^-1
    Biovars.p0v = Biovars.p0./560;% mol N m^-3 s^-1 
end

% Temperature factors
if BioIn.isnem
    Biovars.Kgra = zeros(nbsv,1);
    Biovars.Kgra([Biovars.idx.zs,Biovars.idx.zl,Biovars.idx.zp]) = Np.Kgra(3:5);    
else    
   Biovars.Kgra = zeros(nbsv,1);
   Biovars.Kgra(Biovars.iszoo & ~Biovars.isextrazoo) = Np.Kgra(3:5);
   Biovars.Kgra(Biovars.isextrazoo) = BioIn.kgra;
end

% lambda factor (i.e. ivlev constant) in case I don't use the planktonic
% foraging arena 

if BioIn.isnem
    Biovars.lambda = Biovars.lambda;
else
    Biovars.lambda(Biovars.isextrazoo) = BioIn.lambda;
end


tempfac0 = exp(Biovars.Kgra .* BioIn.temp);
Biovars.gourmetPL = exp(-Biovars.grpusai(2,5).*(BioIn.ZL + BioIn.ZS));
Biovars.gourmetZS = exp(-Biovars.grpusai(3,5).*BioIn.ZL);

Biovars.q0vat0 = bsxfun(@rdivide, Biovars.q0v, tempfac0'); %nb x nb / 1 x nbsv

% Mortality (using ecopath-based mortality)
if BioIn.isnem
  Biovars.m0 = zeros(nbsv,1);
  Biovars.m0(isEwEliving) = Ep.otherMortRate(EwElivingidx); %s-1
else
Biovars.m0 = zeros(nbsv,1);
Biovars.m0(1:EM.nlive) = Ep.otherMortRate(1:EM.nlive);  % s^-1
end

if any(Biovars.m0 < 0)
    warning('BIOMODEL:m0neg', 'Some m0 < 0, setting to 0');
    Biovars.m0 = max(Biovars.m0, 0);
end

% Mortality: Allows for any exponential function aB.^y
if BioIn.isnem
    m0 = zeros(nbsv,1);
    m0(isEwEliving) = Ep.otherMortRate(EwElivingidx); % s-1
else
    m0 = zeros(nbsv,1);
    m0(1:EM.nlive) = Ep.otherMortRate(1:EM.nlive);  % s^-1
end


m0b = m0 .* Biovars.b0; % mass-balanced flux, molN/m^2/s. nbsv x 1 array
                         
                         
Biovars.m0coef = m0b./(Biovars.b0.^Biovars.m0exp);
Biovars.m0coef = max(Biovars.m0coef, 0); % get rid of NaNs (and I guess negatives, if they ever appeared, though they shouldn't)

if Biovars.diapause
    EwEzlidx = strcmp(BioIn.types, 'zl');
    
    bzl = Biovars.b0v(Biovars.idx.zl); %molN/m3
    if BioIn.isnem
        m0bzl = Ep.otherMortRate(EwEzlidx) .* bzl; %molN/m3/s
    else
        m0bzl = Ep.otherMortRate(Biovars.idx.zl) .* bzl; %molN/m3/s
    end
    mexp = Biovars.m0exp(Biovars.idx.zl);
    Biovars.m0coef(Biovars.idx.zl) = 0;
    Biovars.m0coef([Biovars.idx.zl1 Biovars.idx.zl2]) = m0bzl./(bzl.^mexp);
    
    Biovars = setdiapauseparams(BioIn, Biovars, Grd);
end


%**************************************************************************

function [newbio, db, Flx, Diag] = sourcesink(nemflag, oldbio, P, B, G, O2, Arch, it, varargin)

%--------------------------------------------------------------------------
% Set up parameters and forcing data at current time step that will be 
% passed to main ODE function
%--------------------------------------------------------------------------

B.O2 = O2; % nz x nx array current o2


% forcing data at current time step for all relevant critters
% 3 x n cell array. row 1 = pool index; row 2 = time-series type (same code as in Ecosim);
% row 3 = data at current time step; 
if (length(varargin) == 1)
    B.fish = varargin{1};
end

%--------------------------------------------------------------------------
% Split/combine ZL groups as necessary for diapause
%--------------------------------------------------------------------------
if B.diapause
    it = find(B.t == G.t);
    zltot = sum(oldbio(:,:,[B.idx.zl1 B.idx.zl2]), 3);
    if B.zlsplit(it)
        ztransfer = oldbio(:,:, B.idx.zl1) .* B.zlsplit(it);
        oldbio(:,:,B.idx.zl2) = oldbio(:,:,B.idx.zl2) + ztransfer;
        oldbio(:,:,B.idx.zl1) = oldbio(:,:,B.idx.zl1) - ztransfer;
    elseif B.zlcombine(it)
        oldbio(:,:,B.idx.zl1) = zltot;
        oldbio(:,:,B.idx.zl2) = 0;
    end
    oldbio(:,:,B.idx.zl) = 0;
end

%--------------------------------------------------------------------------
% Integrate biology over this time step
%--------------------------------------------------------------------------
[newbio, db, Flx, Diag, badthings] = integratebio(@bioode, nemflag, G, ...
    oldbio, P, B, Arch, it, B.odesolver{:});

if B.diapause
    newbio(:,:,B.idx.zl) = newbio(:,:,B.idx.zl1) + newbio(:,:,B.idx.zl2);
end

% Check and correct for silica issue

[nz,nx,nb] = size(newbio);
A = false(nz,nx,nb);
isneg = newbio(:,:,B.idx.plsi) < 0;

A(:,:,B.idx.plsi) = isneg;

Diag.extrasi = zeros(size(isneg));
C = newbio(:,:,B.idx.plsi);
Diag.extrasi(isneg) = -C(isneg);
newbio(A) = 0;
badthings(A) = false;

% Check for problems

%if any(badthings(:))
%    [ridx,cidx,tidx] = ind2sub(size(badthings),find(badthings));
%    nb = length(ridx);
%    errstr = cell(nb,1);
%    for ii = 1:nb
%        badbio = newbio(ridx(ii), cidx(ii), tidx(ii));
%        if isnan(badbio)
%            errstr{ii} = sprintf('NaN: depth %d, longitude %d, critter %d, time = %d', ridx(ii), cidx(ii), tidx(ii), G.t);
%        elseif isinf(badbio)
%            errstr{ii} = sprintf('Inf: depth %d, longitude %d, critter %d, time = %d', ridx(ii), cidx(ii), tidx(ii), G.t);
%        elseif badbio < 0
%            errstr{ii} = sprintf('Neg: depth %d, longitude %d, critter %d, time = %d', ridx(ii), cidx(ii), tidx(ii), G.t); 
%        end
%    end
%    errstr = sprintf('  %s\n', errstr{:});
%    errstr = sprintf('Biology out of range:\n%s', errstr);
%    error('WCE:biologyOutOfRange', errstr);
%end
                      
%**************************************************************************

function wsink = vertmove(B, G)

wsink = B.settle;

if B.diapause
    
    it = find(B.t == G.t);
    
    if B.zlswim(it) == 1
        wsink(2:end, 2, B.idx.zl2) = 80./86400; % swim up
    elseif B.zlswim(it) == -1
        wsink(1:end-1,2, B.idx.zl2) = -80./86400; % stay down
    end
    
end

