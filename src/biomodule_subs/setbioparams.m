function [Biovars, Np] = setbioparams(BioIn, nemidx, nbsv, Grd, Biovars)
%SETBIOPARAMS Assign and transform NEMURO parameters for the biomodel module
%
% [Biovars, Np] = setbioparams(BioIn, nemidx, nbsv, Grd, Biovars)
%
% This function extracts NEMURO parameters from the input structure, 
% converts units to SI when necessary, and maps them to the correct 
% state variable indices in the Biovars structure for use in the 
% mixed-layer ecological model.
%
% Inputs:
%   BioIn   : Struct containing biological and model configuration, including .NemParam
%   nemidx  : Indices of NEMURO variables in the full state variable list
%   nbsv    : Total number of biological state variables
%   Grd     : Grid structure with spatial information (nx, nz)
%   Biovars : Existing biological variable structure to update
%
% Outputs:
%   Biovars : Updated structure containing all derived biological parameters
%   Np      : Struct of NEMURO parameters reorganized into arrays (from nemparams2arrays)
%
% This file was derived from the original setwcenemparams routine
% developed by Kelly Kearney for the WCE/NEMURO framework and modified
% for the WCVI-E2E coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2008–2015 Kelly Kearney
%
% Modifications by Virginie Bornarel (2017–2026) include:
%   - revised parameter mapping for the WCVI-E2E state-variable structure
%   - addition of oxygen, nitrification, and denitrification parameters
%   - adaptation of settling velocities to a 2D model grid
%   - removal of iron/scavenging and prey-visibility parameters not used
%     in the WCVI-E2E configuration
%   - updated unit conversions and documentation
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

%----------------------------------
% 1. ODE Solver Method
%----------------------------------

if ischar(BioIn.odesolver)
    Biovars.odesolver = {BioIn.odesolver};
else
    Biovars.odesolver = BioIn.odesolver;
end

%----------------------------------
% 2. Get Structured Parameters
%----------------------------------  

Np = nemparams2arrays(BioIn.NemParam);

%----------------------------------
% 3. Scalar Assignments
%----------------------------------

Biovars.alpha1                  = Np.alpha1;              % m^-1  
Biovars.alpha2                  = Np.alpha2/1000;         % (m^3/molN)/m
Biovars.usesteele               = Np.usesteele;           % logical
Biovars.RSiN                    = Np.RSiN;                % molSi/molN
Biovars.RNo3PON                 = Np.RNo3PON;             % molNo3/molPON
Biovars.I0                      = Np.I0 / 0.001433;       % W.m^-2
Biovars.KI                      = Np.KI / 0.001433;       % W.m^-2
Biovars.alphao2                 = Np.alphao2 / 1000;      % m3.mol-1
Biovars.alphaNo3                = Np.alphaNo3 / 1000;     % m3.mol-1
Biovars.alphaNH4                = Np.alphaNH4 / 1000;     % m3.mol-1
Biovars.Nit0                    = Np.Nit0;                % s^-1
Biovars.KNit                    = Np.Knit;                % deg C^-1

%----------------------------------
% 4. Vector Assignments (to nbsv-length arrays)
%----------------------------------
Biovars.Iopt                    = zeros(nbsv,1);
Biovars.Iopt(nemidx(1:11))      = Np.Iopt / 0.001433;     % W/m^2   
Biovars.mor0                    = zeros(nbsv,1);
Biovars.mor0(nemidx(1:11))      = Np.mor0 .*0.001;        % m3/molN/s
Biovars.Kmor                    = zeros(nbsv,1);
Biovars.Kmor(nemidx(1:11))      = Np.Kmor;                % /degC
Biovars.alpha                   = zeros(nbsv,1);
Biovars.alpha(nemidx(1:11))     = Np.alpha * 0.001433;    % (W/m^2)^-1
Biovars.Vmax                    = zeros(nbsv,1);
Biovars.Vmax(nemidx(1:11))      = Np.Vmax;                % s^-1
Biovars.Kno3                    = zeros(nbsv,1);
Biovars.Kno3(nemidx(1:11))      = Np.Kno3 * 1000;         % molN/m^3
Biovars.pusai                   = zeros(nbsv,1);
Biovars.pusai(nemidx(1:11))     = Np.pusai / 1000;        % m^3/molN
Biovars.Knh4                    = zeros(nbsv,1);
Biovars.Knh4(nemidx(1:11))      = Np.Knh4 * 1000;         % molN/m^3 
Biovars.Ksi                     = zeros(nbsv,1);
Biovars.Ksi(nemidx(1:11))       = Np.Ksi * 1000;          % molSi/m^3
Biovars.alphaeg                 = zeros(nbsv,1);
Biovars.alphaeg(nemidx(1:11))   = Np.alphaeg;             % no unit
Biovars.beta                    = zeros(nbsv,1);
Biovars.beta(nemidx(1:11))      = Np.beta;                % no unit
Biovars.Kgpp                    = zeros(nbsv,1);
Biovars.Kgpp(nemidx(1:11))      = Np.Kgpp;                % degC^-1
Biovars.gamma                   = zeros(nbsv,1);
Biovars.gamma(nemidx(1:11))     = Np.gamma;               % no unit
Biovars.lambda                  = zeros(nbsv,1);
Biovars.lambda(nemidx(1:11))    = Np.lambda / 1000;       % m^3/molN
Biovars.res0                    = zeros(nbsv,1);
Biovars.res0(nemidx(1:11))      = Np.res0;                % s^-1
Biovars.Kres                    = zeros(nbsv,1);
Biovars.Kres(nemidx(1:11))      = Np.Kres;                % degC^-1

%----------------------------------
% 5. Matrix Assignments (nbsv x nbsv)
%----------------------------------
Biovars.grmax                   = zeros(nbsv,nbsv);
Biovars.grmax(nemidx(1:11),nemidx(1:11)) = Np.grmax;      % s^-1
Biovars.thresh                  = zeros(nbsv,nbsv);
Biovars.thresh(nemidx(1:11),nemidx(1:11)) = Np.thresh * 1000; % molN/m^3
Biovars.grpusai                 = zeros(nbsv,nbsv);
Biovars.grpusai(nemidx(1:11),nemidx(1:11)) = Np.grpusai / 1000; %m^3/molN
Biovars.Kdec                    = zeros(nbsv,nbsv);
Biovars.Kdec(nemidx(1:11),nemidx(1:11)) = Np.Kdec;        % degC^-1
Biovars.vdec                    = zeros(nbsv,nbsv);
Biovars.vdec(nemidx(1:11),nemidx(1:11)) = Np.vdec;        % s^-1

%----------------------------------
% 6. Settling Velocities (optional sink term)
%----------------------------------
Biovars.settle = zeros(Grd.nz, Grd.nx, nbsv);
Biovars.settle(:,:,nemidx(1:11)) = repmat(reshape(-Np.settle,1,1,[]), Grd.nz, Grd.nx);  % m/s,

