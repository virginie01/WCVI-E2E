function In = loadPhysicalForcing(Forcingdir, LBdir)
%LOADPHYSICALFORCING Load physical forcing datasets from .mat files
%   Forcingdir : directory containing most forcing files
%   LBdir      : directory containing LBtmp_input and LBsal_input (lateral
%                boundary conditions)

forcingFiles = {
    'Qi_input'
    'airtmp_input'
    'dewptT_input'
    'uWndSpd_input'
    'vWndSpd_input'
    'xfilphy_input'
    'xfil_input'
    'tauy2_input'
    'mld_input'
    'entrnmnt_input'
    'p_input'
    'cu_input'
    'dc_input'
    'sbc_input'
};

lbFiles = {
    'LBtmp_input'
    'LBsal_input'
};

% Load regular forcing files
for i = 1:numel(forcingFiles)
    fname = forcingFiles{i};
    S = load(fullfile(Forcingdir, [fname '.mat']));
    In.(fname) = S.(fname);
end

% Load LB files
for i = 1:numel(lbFiles)
    fname = lbFiles{i};
    S = load(fullfile(LBdir, [fname '.mat']));
    In.(fname) = S.(fname);
end

end