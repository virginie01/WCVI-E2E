function Biovars = setdiapauseparams(BioIn, Biovars, Grd)
%SETDIAPAUSEPARAMS Sets seasonal migration (diapause) parameters for large zooplankton (ZL)
%
% Biovars = setdiapauseparams(BioIn, Biovars, Grd)
%
% This function configures seasonal vertical migration of large zooplankton (ZL),
% splitting them into two subgroups (ZL1 active and ZL2 dormant) and tracking their 
% behavior over time according to model inputs.
%
% Outputs:
%   Biovars - Updated struct with:
%
%             .zlswim:   [nt+1 x 1]. Indicates 1, 0 or -1 depending on 
%                        whether the current time step corresponds to upward 
%                        migration (i.e. between Sep 15 and April 1), no 
%                        directed movement (i.e. between April 1 and May 1 )
%                        or downward migration (i.e. between May 1 and Sep 15)
%                        Note: no directed movement before the first downward migration
%                        From January 1992 to Sept 1992
%             .zlsplit:  [nt+1 x 1]. Percent of ZL1 to transfer to ZL2 at 
%                        each time step during the whole simulation (i.e.
%                        7% on shelf and 12% on slope between May 1 and Sep 1 and 0% otherwise)
%
%             .zlcombine:[nt x 1] logical, indicate at what time step
%                        both ZL2 and ZL1 recombine (i.e. April 1)
%
%             .t:        [1 x nt]. Equivalent to Grd.time, time elapsed 
%                        from model start time to the beginning of each time 
%                        interval (seconds) 
%
% This file was derived from the original setdiapauseparams routine
% developed by Kelly Kearney for the WCE/NEMURO framework and modified
% for the WCVI-E2E coastal upwelling ecosystem model.
%
% Original framework:
% Copyright (c) 2014 Kelly Kearney
%
% Modifications and extensions by Virginie Bornarel (2017–2026) include:
%   - revised diapause parameter handling for WCVI-E2E biological groups
%   - updated seasonal migration timing and split/recombination logic
%   - support for WCVI-E2E-specific diapause parameterization
%   - expanded documentation and ecological interpretation of ZL states
%   - revised year-handling logic for seasonal transitions
%
% Distributed under the MIT License.
% See LICENSE file in the repository root for details.

%% Transfer ZL parameters to ZL1 and ZL2
if BioIn.isnem
    epvars = {'alphaeg', 'beta', 'mor0', 'Kmor'};
else
    epvars = {'gs', 'ge', 'm0exp'};
end
for ii = 1:length(epvars)
    Biovars.(epvars{ii})([Biovars.idx.zl1 Biovars.idx.zl2]) = Biovars.(epvars{ii})(Biovars.idx.zl);
end

%% Setup migration timing

yr = Grd.start_date(1):Grd.end_date(1);
yr = [yr yr(end)+1];

d1 = datevec(BioIn.dpEnd);
d2 = datevec(datenum(d1) + BioIn.dpEndSpan);
d3 = datevec(BioIn.dpStart);

bin = [];
val = [];
for iy = 1:length(yr)-1
    bintmp = [d1; d2; d3];
    if d1(1) == 2021 && d2(1) == 2021
    bintmp(:,1) = [yr(iy);yr(iy);yr(iy+1)];
    elseif d1(1)==2022 && d2(1)==2022
    bintmp(:,1) = [yr(iy);yr(iy);yr(iy)];
    elseif d1(1)==2021 && d2(1)==2022
    bintmp(:,1) = [yr(iy);yr(iy+1);yr(iy+1)];
    end        
    bin = [bin; bintmp];
    val = [val; 1; 0; -1];
end
if datenum(Grd.start_date) < datenum(bin(1,:))
    bin = [Grd.start_date; bin];
    val = [0; val];
end
if datenum(Grd.end_date) > datenum(bin(end,:))
    bin = [bin; Grd.end_date];
    val = [val; val(end)];
end
bin = datenum(bin);

dnsim = datenum(Grd.start_date) + Grd.time/86400;
[~, binidx] = histc(dnsim, bin);


Biovars.zlswim = val(binidx);
idx = find(Biovars.zlswim < 0, 1);
Biovars.zlswim(1:idx-1) = 0; 


%% Setup ZL1 → ZL2 split fraction
d1 = datevec(BioIn.dpStart);
d2 = datevec(datenum(d1) + BioIn.dpSpan);

frac = BioIn.dpPercent./100;

bin = [];
val = [];
for iy = 1:length(yr)
    bintmp = [d1; d2];
    bintmp(:,1) = yr(iy);
    bin = [bin; bintmp];
    val = [val; frac; 0];
end
if datenum(Grd.start_date) < datenum(bin(1,:))
    bin = [Grd.start_date; bin];
    val = [0; val];
end
if datenum(Grd.end_date) > datenum(bin(end,:))
    bin = [bin; Grd.end_date];
    val = [val; val(end)];
end
bin = datenum(bin);
 
dnsim = datenum(Grd.start_date) + Grd.time/86400;
[~, binidx] = histc(dnsim, bin);
Biovars.zlsplit = val(binidx);

%% Setup ZL1/ZL2 recombination marker
Biovars.zlcombine = [false; Biovars.zlswim(2:end) == 0 & Biovars.zlswim(1:end-1) == 1];

%% Pass time vector for convenience
Biovars.t = Grd.time;


