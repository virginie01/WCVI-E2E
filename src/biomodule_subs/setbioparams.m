function [Biovars, Np] = setbioparams(BioIn, nemidx, nbsv, Grd, Biovars)
%SETWCENEMPARAMS Setup of parameters for mixed_layer wce/nemurokak modules

%-----------------------
% Shared parameters
%-----------------------

% ODE solver

if ischar(BioIn.odesolver)
    Biovars.odesolver = {BioIn.odesolver};
else
    Biovars.odesolver = BioIn.odesolver;
end

% Extended NEMURO parameters with 1:1 correspondence to nemurokak
% (only defined for original-11 NEMURO-derived variables, so reassign to
% proper indices) 

Np = nemparams2arrays(BioIn.NemParam);

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

% PON and Opal sink

% if BioIn.isnem
%     Biovars.sink = zeros(1,nbsv);
%     Biovars.sink(1,Biovars.idx.pon)  = -Np.settle(Biovars.idx.pon); 
%     Biovars.sink(1,Biovars.idx.opal) = -Np.settle(Biovars.idx.opal);
%     Biovars.sink(1,Biovars.idx.pofe) = -Np.settle(Biovars.idx.pon); % Same as PON
% else
Biovars.settle = zeros(Grd.nz, Grd.nx, nbsv);
Biovars.settle(:,:,nemidx(1:11)) = repmat(reshape(-Np.settle,1,1,[]), Grd.nz, Grd.nx);  % m/s,

