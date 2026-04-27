function BIOTS = buildFisheriesTS(fisheriesTS_input, Grd)
%BUILDFISHERIESTS Construct Ecosim-style fisheries time series structure
%
% Output matches legacy BIOTS structure exactly.

% Extract time rows
dv = cell2mat(fisheriesTS_input(5:end,1:6));
dn = datenum(dv);

dnstart = datenum(Grd.start_date);
tsec = (dn - dnstart) * 86400;

% Convert to cell array (legacy behavior)
t = num2cell(tsec);

% Header rows (Title / Weight / Pool / Type)
col = fisheriesTS_input(3:4,7:end);

% Data matrix
data = fisheriesTS_input(5:end,7:end);

BIOTS = struct( ...
    't',    {t}, ...
    'col',  {col}, ...
    'data', {data} );

end