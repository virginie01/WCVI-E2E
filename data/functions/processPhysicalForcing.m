function Forcing = processPhysicalForcing(In, Grd)
%PROCESSPHYSICALFORCING Format, interpolate, smooth, and derive forcings

%% Formatting: 

Qi      = initinterpdata('time and surface grid', In.Qi_input, Grd);
airtmp  = initinterpdata('time and surface grid', In.airtmp_input, Grd);
dewptT  = initinterpdata('time and surface grid', In.dewptT_input, Grd);

UWndSpd = initinterpdata('time and surface grid', In.uWndSpd_input, Grd);
VWndSpd = initinterpdata('time and surface grid', In.vWndSpd_input, Grd);

Xfil    = initinterpdata('time', In.xfil_input, Grd);
Xfilphy = initinterpdata('time', In.xfilphy_input, Grd);
Tauy2   = initinterpdata('time', In.tauy2_input, Grd);

Mld     = initinterpdata('time and surface grid', In.mld_input, Grd);
Entrnmnt= initinterpdata('time and surface grid', In.entrnmnt_input, Grd);

P       = initinterpdata('time and surface grid', In.p_input, Grd);
CU      = initinterpdata('time', In.cu_input, Grd);
DC      = initinterpdata('time', In.dc_input, Grd);
SBC     = initinterpdata('time', In.sbc_input, Grd);

tmp     = initinterpdata('time and boundaries', In.LBtmp_input, Grd);
sal     = initinterpdata('time and boundaries', In.LBsal_input, Grd);

%% Interpolation to the ODE solver timestep

interpToDataTime = @(A) struct( ...
    't', Grd.datatime, ...
    'data', interp1(A.t, A.data, Grd.datatime) ...
);

Qi      = interpToDataTime(Qi);
airtmp  = interpToDataTime(airtmp);
dewptT  = interpToDataTime(dewptT);
UWndSpd = interpToDataTime(UWndSpd);
VWndSpd = interpToDataTime(VWndSpd);
Xfil    = interpToDataTime(Xfil);
Xfilphy = interpToDataTime(Xfilphy);
Tauy2   = interpToDataTime(Tauy2);
Mld     = interpToDataTime(Mld);
Entrnmnt= interpToDataTime(Entrnmnt);
P       = interpToDataTime(P);
CU      = interpToDataTime(CU);
DC      = interpToDataTime(DC);
SBC     = interpToDataTime(SBC);
tmp     = interpToDataTime(tmp);
sal     = interpToDataTime(sal);

%% Smoothing - daily mean for time series with diurnal fluctuations

span = (24.*3600) ./ In.datadt;

for i = 1:Grd.nx
    Qi.data(:,i)     = smooth(Qi.data(:,i), span);
    airtmp.data(:,i)= smooth(airtmp.data(:,i), span);
    dewptT.data(:,i)= smooth(dewptT.data(:,i), span);
    UWndSpd.data(:,i)= smooth(UWndSpd.data(:,i), span);
    VWndSpd.data(:,i)= smooth(VWndSpd.data(:,i), span);
    P.data(:,i)      = smooth(P.data(:,i), span);
end

%% Derived physics

Qo.data = clearsky(Grd.start_date, Grd.datatime, In.Lat)';
Qo.t = Grd.datatime;

[~, tauy, u10, v10] = wstress(UWndSpd.data, VWndSpd.data, In.whgt);
wspeed10 = abs(u10 + sqrt(-1)*v10);

% Translate from dynes/cm2 to Newtons/m2

tauy = 0.1*tauy;

Tauy.t = UWndSpd.t;
Tauy.x = UWndSpd.x;
Tauy.data = tauy;

Wspd10.t = UWndSpd.t;
Wspd10.x = UWndSpd.x;
Wspd10.data = wspeed10;

%% Package outputs

t = Xfilphy.t;
dnstart = datenum(Grd.start_date);
dncurrent = dnstart + (t./86400);
t = datevec(dncurrent);

QI = [t Qi.data];
AIRTMP = [t airtmp.data];
DEWPTT = [t dewptT.data];
QO = [t Qo.data];
TAUY = [t Tauy.data];
TAUY2 = [t Tauy2.data];
WSPD10 = [t Wspd10.data];
XFIL = [t Xfil.data];
XFILPHY = [t Xfilphy.data];
MLD = [t Mld.data];
ENT = [t Entrnmnt.data];
TMP = [t tmp.data];
SAL = [t sal.data];
PRATE = [t P.data];
CU = [t CU.data];
DC = [t DC.data];
SBC = [t SBC.data];

Forcing = struct();
Forcing.QI = QI;
Forcing.AIRTMP = AIRTMP;
Forcing.DEWPTT = DEWPTT;
Forcing.QO = QO;
Forcing.TAUY = TAUY;
Forcing.TAUY2 = TAUY2;
Forcing.WSPD10 = WSPD10;
Forcing.XFIL = XFIL;
Forcing.XFILPHY = XFILPHY;
Forcing.MLD = MLD;
Forcing.ENT = ENT;
Forcing.TMP = TMP;
Forcing.SAL = SAL;
Forcing.PRATE = PRATE;
Forcing.CU = CU;
Forcing.DC = DC;
Forcing.SBC = SBC;
end

