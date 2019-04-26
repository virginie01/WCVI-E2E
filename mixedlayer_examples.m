%% Examples for the mixed_layer model

%% Setting up your path
% 
% The mixed_layer code is spread over several folders (mostly because that
% makes it easier to maintain on my end... sorry about that).  To get the
% basics running, you need to add a few folders to your path:  
% 
addpath('D:/These UBC/WCVI-E2E/physical_model');
addpath('D:/These UBC/WCVI-E2E/mergestruct');
addpath('D:/These UBC/WCVI-E2E/seawater_ver3_2');
%
% (Alter the paths as necessary to wherever you place these folders).

%% Running physics-only simulations
%
% The simplest way to run the |WCVIE2E_physicalmodel| model is to do so without any
% of the biological modules turned on.  The code comes with some default forcing datasets which can be
% used for some quick tests.  To run with all default variables, simply
% specify an output folder name.  

WCVIE2E_physicalmodel('default');

% When run in single-simulation mode like this, the model creates two
% output files: one that holds the dimension variables (dimensions.nc), and
% one for everything else (sim0001.nc).  

% add new path: folder containing output files 

addpath('D:/These UBC/WCVI-E2E/default');

% read in output variables

Tmp=ncread('sim0001.nc','temp');
Sal=ncread('sim0001.nc','sal');
Sig=ncread('sim0001.nc','sig');

%To run with diagnostic tests:

WCVIE2E_physicalmodel('default','Tdiag',true,'Sdiag',true);

% add new path: folder containing diagnostic output files 

addpath('D:/These UBC/WCVI-E2E/T_diagnostic');
addpath('D:/These UBC/WCVI-E2E/Sal_diagnostic');

%%%%%%%%%%---------------------TEMPERATURE--------------------%%%%%%%%%%%%%
% read in output variables
tmp=ncread(fullfile('T_diagnostic','data_values.nc'),'temp');
sig=ncread(fullfile('T_diagnostic','data_values.nc'),'sig');
k1T=ncread(fullfile('T_diagnostic','data_values.nc'),'ode1T');
V1T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeV1T');
H1T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeH1T');
X1T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeX1T');
k2T=ncread(fullfile('T_diagnostic','data_values.nc'),'ode2T');
V2T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeV2T');
H2T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeH2T');
X2T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeX2T');
k3T=ncread(fullfile('T_diagnostic','data_values.nc'),'ode3T');
V3T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeV3T');
H3T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeH3T');
X3T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeX3T');
k4T=ncread(fullfile('T_diagnostic','data_values.nc'),'ode4T');
V4T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeV4T');
H4T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeH4T');
X4T=ncread(fullfile('T_diagnostic','data_values.nc'),'odeX4T');

%--------------------plotting diagnostic outputs-------------------------%

%%%-----temperature-----%
tmp11=tmp(1,1,:);
tmp11=tmp11(:);

tmp21=tmp(2,1,:);
tmp21=tmp21(:);

tmp31=tmp(3,1,:);
tmp31=tmp31(:);

tmp12=tmp(1,2,:);
tmp12=tmp12(:);

tmp22=tmp(2,2,:);
tmp22=tmp22(:);

tmp32=tmp(3,2,:);
tmp32=tmp32(:);

subplot(3,2,1)
plot(tmp11)
title('shelf UL')
xlabel('time step')
ylabel('temperature deg')

subplot(3,2,2)
plot(tmp12)
title('slope UL')
xlabel('time step')
ylabel('temperature deg')

subplot(3,2,3)
plot(tmp21)
title('shelf LL')
xlabel('time step')
ylabel('temperature deg')

subplot(3,2,4)
plot(tmp22)
title('slope LL')
xlabel('time step')
ylabel('temperature deg')

subplot(3,2,5)
plot(tmp31)
title('shelf Dem')
xlabel('time step')
ylabel('temperature deg')

subplot(3,2,6)
plot(tmp32)
title('slope Dem')
xlabel('time step')
ylabel('temperature deg')

%%%-----density-----%
sig11=sig(1,1,:);
sig11=sig11(:);

sig21=sig(2,1,:);
sig21=sig21(:);

sig31=sig(3,1,:);
sig31=sig31(:);

sig12=sig(1,2,:);
sig12=sig12(:);

sig22=sig(2,2,:);
sig22=sig22(:);

sig32=sig(3,2,:);
sig32=sig32(:);

subplot(3,2,1)
plot(sig11)
title('shelf UL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,2)
plot(sig12)
title('slope UL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,3)
plot(sig21)
title('shelf LL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,4)
plot(sig22)
title('slope LL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,5)
plot(sig31)
title('shelf Dem')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,6)
plot(sig32)
title('slope Dem')
xlabel('time step')
ylabel('density kg.m-3')

%%%-----k1T,k2T,k3T,k4T-----%

k1T11=k1T(1,1,:);k2T11=k2T(1,1,:);k3T11=k3T(1,1,:);k4T11=k4T(1,1,:);
k1T11=k1T11(:);k2T11=k2T11(:);k3T11=k3T11(:);k4T11=k4T11(:);

k1T21=k1T(2,1,:);k2T21=k2T(2,1,:);k3T21=k3T(2,1,:);k4T21=k4T(2,1,:);
k1T21=k1T21(:);k2T21=k2T21(:);k3T21=k3T21(:);k4T21=k4T21(:);

k1T31=k1T(3,1,:);k2T31=k2T(3,1,:);k3T31=k3T(3,1,:);k4T31=k4T(3,1,:);
k1T31=k1T31(:);k2T31=k2T31(:);k3T31=k3T31(:);k4T31=k4T31(:);

k1T12=k1T(1,2,:);k2T12=k2T(1,2,:);k3T12=k3T(1,2,:);k4T12=k4T(1,2,:);
k1T12=k1T12(:);k2T12=k2T12(:);k3T12=k3T12(:);k4T12=k4T12(:);

k1T22=k1T(2,2,:);k2T22=k2T(2,2,:);k3T22=k3T(2,2,:);k4T22=k4T(2,2,:);
k1T22=k1T22(:);k2T22=k2T22(:);k3T22=k3T22(:);k4T22=k4T22(:);

k1T32=k1T(3,2,:);k2T32=k2T(3,2,:);k3T32=k3T(3,2,:);k4T32=k4T(3,2,:);
k1T32=k1T32(:);k2T32=k2T32(:);k3T32=k3T32(:);k4T32=k4T32(:);

subplot(3,2,1)
plot(k1T11)
title('shelf UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(k2T11);plot(k3T11);plot(k4T11);
line([0, 40912],[0,0], 'Color', 'k');
legend('k1T','k2T','k3T','k4T')
hold off

subplot(3,2,2)
plot(k1T12)
title('slope UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(k2T12); plot(k3T12); plot(k4T12)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1T','k2T','k3T','k4T')
hold off

subplot(3,2,3)
plot(k1T21)
title('shelf LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(k2T21); plot(k3T21); plot(k4T21)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1T','k2T','k3T','k4T')
hold off

subplot(3,2,4)
plot(k1T22)
title('slope LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(k2T22); plot(k3T22); plot(k4T22)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1T','k2T','k3T','k4T')
hold off

subplot(3,2,5)
plot(k1T31)
title('shelf Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(k2T31); plot(k3T31); plot(k4T31)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1T','k2T','k3T','k4T')
hold off

subplot(3,2,6)
plot(k1T32)
title('slope Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(k2T32); plot(k3T32); plot(k4T32)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1T','k2T','k3T','k4T')
hold off

%%%-----V1T,V2T,V3T,V4T------%

V1T11=V1T(1,1,:);V2T11=V2T(1,1,:);V3T11=V3T(1,1,:);V4T11=V4T(1,1,:);
V1T11=V1T11(:);V2T11=V2T11(:);V3T11=V3T11(:);V4T11=V4T11(:);

V1T21=V1T(2,1,:);V2T21=V2T(2,1,:);V3T21=V3T(2,1,:);V4T21=V4T(2,1,:);
V1T21=V1T21(:);V2T21=V2T21(:);V3T21=V3T21(:);V4T21=V4T21(:);

V1T31=V1T(3,1,:);V2T31=V2T(3,1,:);V3T31=V3T(3,1,:);V4T31=V4T(3,1,:);
V1T31=V1T31(:);V2T31=V2T31(:);V3T31=V3T31(:);V4T31=V4T31(:);

V1T12=V1T(1,2,:);V2T12=V2T(1,2,:);V3T12=V3T(1,2,:);V4T12=V4T(1,2,:);
V1T12=V1T12(:);V2T12=V2T12(:);V3T12=V3T12(:);V4T12=V4T12(:);

V1T22=V1T(2,2,:);V2T22=V2T(2,2,:);V3T22=V3T(2,2,:);V4T22=V4T(2,2,:);
V1T22=V1T22(:);V2T22=V2T22(:);V3T22=V3T22(:);V4T22=V4T22(:);

V1T32=V1T(3,2,:);V2T32=V2T(3,2,:);V3T32=V3T(3,2,:);V4T32=V4T(3,2,:);
V1T32=V1T32(:);V2T32=V2T32(:);V3T32=V3T32(:);V4T32=V4T32(:);

subplot(3,2,1)
plot(V1T11)
title('shelf UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(V2T11);plot(V3T11);plot(V4T11);
line([0, 40912],[0,0], 'Color', 'k');
legend('V1T','V2T','V3T','V4T')
hold off

subplot(3,2,2)
plot(V1T12)
title('slope UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(V2T12); plot(V3T12); plot(V4T12)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1T','V2T','V3T','V4T')
hold off

subplot(3,2,3)
plot(V1T21)
title('shelf LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(V2T21); plot(V3T21); plot(V4T21)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1T','V2T','V3T','V4T')
hold off

subplot(3,2,4)
plot(V1T22)
title('slope LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(V2T22); plot(V3T22); plot(V4T22)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1T','V2T','V3T','V4T')
hold off

subplot(3,2,5)
plot(V1T31)
title('shelf Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(V2T31); plot(V3T31); plot(V4T31)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1T','V2T','V3T','V4T')
hold off

subplot(3,2,6)
plot(V1T32)
title('slope Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(V2T32); plot(V3T32); plot(V4T32)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1T','V2T','V3T','V4T')
hold off

%%%-----H1T,H2T,H3T,H4T------%

H1T11=H1T(1,1,:);H2T11=H2T(1,1,:);H3T11=H3T(1,1,:);H4T11=H4T(1,1,:);
H1T11=H1T11(:);H2T11=H2T11(:);H3T11=H3T11(:);H4T11=H4T11(:);

H1T21=H1T(2,1,:);H2T21=H2T(2,1,:);H3T21=H3T(2,1,:);H4T21=H4T(2,1,:);
H1T21=H1T21(:);H2T21=H2T21(:);H3T21=H3T21(:);H4T21=H4T21(:);

H1T31=H1T(3,1,:);H2T31=H2T(3,1,:);H3T31=H3T(3,1,:);H4T31=H4T(3,1,:);
H1T31=H1T31(:);H2T31=H2T31(:);H3T31=H3T31(:);H4T31=H4T31(:);

H1T12=H1T(1,2,:);H2T12=H2T(1,2,:);H3T12=H3T(1,2,:);H4T12=H4T(1,2,:);
H1T12=H1T12(:);H2T12=H2T12(:);H3T12=H3T12(:);H4T12=H4T12(:);

H1T22=H1T(2,2,:);H2T22=H2T(2,2,:);H3T22=H3T(2,2,:);H4T22=H4T(2,2,:);
H1T22=H1T22(:);H2T22=H2T22(:);H3T22=H3T22(:);H4T22=H4T22(:);

H1T32=H1T(3,2,:);H2T32=H2T(3,2,:);H3T32=H3T(3,2,:);H4T32=H4T(3,2,:);
H1T32=H1T32(:);H2T32=H2T32(:);H3T32=H3T32(:);H4T32=H4T32(:);

subplot(3,2,1)
plot(H1T11)
title('shelf UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(H2T11);plot(H3T11);plot(H4T11);
line([0, 40912],[0,0], 'Color', 'k');
legend('H1T','H2T','H3T','H4T')
hold off

subplot(3,2,2)
plot(H1T12)
title('slope UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(H2T12); plot(H3T12); plot(H4T12)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1T','H2T','H3T','H4T')
hold off

subplot(3,2,3)
plot(H1T21)
title('shelf LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(H2T21); plot(H3T21); plot(H4T21)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1T','H2T','H3T','H4T')
hold off

subplot(3,2,4)
plot(H1T22)
title('slope LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(H2T22); plot(H3T22); plot(H4T22)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1T','H2T','H3T','H4T')
hold off

subplot(3,2,5)
plot(H1T31)
title('shelf Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(H2T31); plot(H3T31); plot(H4T31)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1T','H2T','H3T','H4T')
hold off

subplot(3,2,6)
plot(H1T32)
title('slope Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(H2T32); plot(H3T32); plot(H4T32)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1T','H2T','H3T','H4T')
hold off

%%%-----X1T,X2T,X3T,X4T------%

X1T11=X1T(1,1,:);X2T11=X2T(1,1,:);X3T11=X3T(1,1,:);X4T11=X4T(1,1,:);
X1T11=X1T11(:);X2T11=X2T11(:);X3T11=X3T11(:);X4T11=X4T11(:);

X1T21=X1T(2,1,:);X2T21=X2T(2,1,:);X3T21=X3T(2,1,:);X4T21=X4T(2,1,:);
X1T21=X1T21(:);X2T21=X2T21(:);X3T21=X3T21(:);X4T21=X4T21(:);

X1T31=X1T(3,1,:);X2T31=X2T(3,1,:);X3T31=X3T(3,1,:);X4T31=X4T(3,1,:);
X1T31=X1T31(:);X2T31=X2T31(:);X3T31=X3T31(:);X4T31=X4T31(:);

X1T12=X1T(1,2,:);X2T12=X2T(1,2,:);X3T12=X3T(1,2,:);X4T12=X4T(1,2,:);
X1T12=X1T12(:);X2T12=X2T12(:);X3T12=X3T12(:);X4T12=X4T12(:);

X1T22=X1T(2,2,:);X2T22=X2T(2,2,:);X3T22=X3T(2,2,:);X4T22=X4T(2,2,:);
X1T22=X1T22(:);X2T22=X2T22(:);X3T22=X3T22(:);X4T22=X4T22(:);

X1T32=X1T(3,2,:);X2T32=X2T(3,2,:);X3T32=X3T(3,2,:);X4T32=X4T(3,2,:);
X1T32=X1T32(:);X2T32=X2T32(:);X3T32=X3T32(:);X4T32=X4T32(:);

subplot(3,2,1)
plot(X1T11)
title('shelf UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(X2T11);plot(X3T11);plot(X4T11);
line([0, 40912],[0,0], 'Color', 'k');
legend('X1T','X2T','X3T','X4T')
hold off

subplot(3,2,2)
plot(X1T12)
title('slope UL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(X2T12); plot(X3T12); plot(X4T12)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1T','X2T','X3T','X4T')
hold off

subplot(3,2,3)
plot(X1T21)
title('shelf LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(X2T21); plot(X3T21); plot(X4T21)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1T','X2T','X3T','X4T')
hold off

subplot(3,2,4)
plot(X1T22)
title('slope LL')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(X2T22); plot(X3T22); plot(X4T22)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1T','X2T','X3T','X4T')
hold off

subplot(3,2,5)
plot(X1T31)
title('shelf Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(X2T31); plot(X3T31); plot(X4T31)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1T','X2T','X3T','X4T')
hold off

subplot(3,2,6)
plot(X1T32)
title('slope Dem')
xlabel('time step')
ylabel('deg C.s^-1')
hold on
plot(X2T32); plot(X3T32); plot(X4T32)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1T','X2T','X3T','X4T')
hold off

%%%%%%%%%%---------------------SALINITY--------------------%%%%%%%%%%%%%%%

% read in output variables
sal=ncread(fullfile('Sal_diagnostic','data_values.nc'),'sal');
sig=ncread(fullfile('Sal_diagnostic','data_values.nc'),'sig');
k1S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'ode1S');
V1S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeV1S');
H1S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeH1S');
X1S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeX1S');
k2S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'ode2S');
V2S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeV2S');
H2S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeH2S');
X2S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeX2S');
k3S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'ode3S');
V3S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeV3S');
H3S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeH3S');
X3S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeX3S');
k4S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'ode4S');
V4S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeV4S');
H4S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeH4S');
X4S=ncread(fullfile('Sal_diagnostic','data_values.nc'),'odeX4S');

%--------------------plotting diagnostic outputs-------------------------%

%%%-----salinity-----%
sal11=sal(1,1,:);
sal11=sal11(:);

sal21=sal(2,1,:);
sal21=sal21(:);

sal31=sal(3,1,:);
sal31=sal31(:);

sal12=sal(1,2,:);
sal12=sal12(:);

sal22=sal(2,2,:);
sal22=sal22(:);

sal32=sal(3,2,:);
sal32=sal32(:);

subplot(3,2,1)
plot(sal11)
title('shelf UL')
xlabel('time step')
ylabel('salinity psu')

subplot(3,2,2)
plot(sal12)
title('slope UL')
xlabel('time step')
ylabel('salinity psu')

subplot(3,2,3)
plot(sal21)
title('shelf LL')
xlabel('time step')
ylabel('salinity psu')

subplot(3,2,4)
plot(sal22)
title('slope LL')
xlabel('time step')
ylabel('salinity psu')

subplot(3,2,5)
plot(sal31)
title('shelf Dem')
xlabel('time step')
ylabel('salinity psu')

subplot(3,2,6)
plot(sal32)
title('slope Dem')
xlabel('time step')
ylabel('salinity psu')

%%%-----density-----%

sig11=sig(1,1,:);
sig11=sig11(:);

sig21=sig(2,1,:);
sig21=sig21(:);

sig31=sig(3,1,:);
sig31=sig31(:);

sig12=sig(1,2,:);
sig12=sig12(:);

sig22=sig(2,2,:);
sig22=sig22(:);

sig32=sig(3,2,:);
sig32=sig32(:);

subplot(3,2,1)
plot(sig11)
title('shelf UL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,2)
plot(sig12)
title('slope UL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,3)
plot(sig21)
title('shelf LL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,4)
plot(sig22)
title('slope LL')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,5)
plot(sig31)
title('shelf Dem')
xlabel('time step')
ylabel('density kg.m-3')

subplot(3,2,6)
plot(sig32)
title('slope Dem')
xlabel('time step')
ylabel('density kg.m-3')


%%%-----k1T,k2T,k3T,k4T-----%

k1S11=k1S(1,1,:);k2S11=k2S(1,1,:);k3S11=k3S(1,1,:);k4S11=k4S(1,1,:);
k1S11=k1S11(:);k2S11=k2S11(:);k3S11=k3S11(:);k4S11=k4S11(:);

k1S21=k1S(2,1,:);k2S21=k2S(2,1,:);k3S21=k3S(2,1,:);k4S21=k4S(2,1,:);
k1S21=k1S21(:);k2S21=k2S21(:);k3S21=k3S21(:);k4S21=k4S21(:);

k1S31=k1S(3,1,:);k2S31=k2S(3,1,:);k3S31=k3S(3,1,:);k4S31=k4S(3,1,:);
k1S31=k1S31(:);k2S31=k2S31(:);k3S31=k3S31(:);k4S31=k4S31(:);

k1S12=k1S(1,2,:);k2S12=k2S(1,2,:);k3S12=k3S(1,2,:);k4S12=k4S(1,2,:);
k1S12=k1S12(:);k2S12=k2S12(:);k3S12=k3S12(:);k4S12=k4S12(:);

k1S22=k1S(2,2,:);k2S22=k2S(2,2,:);k3S22=k3S(2,2,:);k4S22=k4S(2,2,:);
k1S22=k1S22(:);k2S22=k2S22(:);k3S22=k3S22(:);k4S22=k4S22(:);

k1S32=k1S(3,2,:);k2S32=k2S(3,2,:);k3S32=k3S(3,2,:);k4S32=k4S(3,2,:);
k1S32=k1S32(:);k2S32=k2S32(:);k3S32=k3S32(:);k4S32=k4S32(:);

subplot(3,2,1)
plot(k1S11)
title('shelf UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(k2S11);plot(k3S11);plot(k4S11);
line([0, 40912],[0,0], 'Color', 'k');
legend('k1S','k2S','k3S','k4S')
hold off

subplot(3,2,2)
plot(k1S12)
title('slope UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(k2S12); plot(k3S12); plot(k4S12)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1S','k2S','k3S','k4S')
hold off

subplot(3,2,3)
plot(k1S21)
title('shelf LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(k2S21); plot(k3S21); plot(k4S21)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1S','k2S','k3S','k4S')
hold off

subplot(3,2,4)
plot(k1S22)
title('slope LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(k2S22); plot(k3S22); plot(k4S22)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1S','k2S','k3S','k4S')
hold off

subplot(3,2,5)
plot(k1S31)
title('shelf Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(k2S31); plot(k3S31); plot(k4S31)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1S','k2S','k3S','k4S')
hold off

subplot(3,2,6)
plot(k1S32)
title('slope Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(k2S32); plot(k3S32); plot(k4S32)
line([0, 40912],[0,0], 'Color', 'k');
legend('k1S','k2S','k3S','k4S')
hold off

%%%-----V1T,V2T,V3T,V4T------%

V1S11=V1S(1,1,:);V2S11=V2S(1,1,:);V3S11=V3S(1,1,:);V4S11=V4S(1,1,:);
V1S11=V1S11(:);V2S11=V2S11(:);V3S11=V3S11(:);V4S11=V4S11(:);

V1S21=V1S(2,1,:);V2S21=V2S(2,1,:);V3S21=V3S(2,1,:);V4S21=V4S(2,1,:);
V1S21=V1S21(:);V2S21=V2S21(:);V3S21=V3S21(:);V4S21=V4S21(:);

V1S31=V1S(3,1,:);V2S31=V2S(3,1,:);V3S31=V3S(3,1,:);V4S31=V4S(3,1,:);
V1S31=V1S31(:);V2S31=V2S31(:);V3S31=V3S31(:);V4S31=V4S31(:);

V1S12=V1S(1,2,:);V2S12=V2S(1,2,:);V3S12=V3S(1,2,:);V4S12=V4S(1,2,:);
V1S12=V1S12(:);V2S12=V2S12(:);V3S12=V3S12(:);V4S12=V4S12(:);

V1S22=V1S(2,2,:);V2S22=V2S(2,2,:);V3S22=V3S(2,2,:);V4S22=V4S(2,2,:);
V1S22=V1S22(:);V2S22=V2S22(:);V3S22=V3S22(:);V4S22=V4S22(:);

V1S32=V1S(3,2,:);V2S32=V2S(3,2,:);V3S32=V3S(3,2,:);V4S32=V4S(3,2,:);
V1S32=V1S32(:);V2S32=V2S32(:);V3S32=V3S32(:);V4S32=V4S32(:);

subplot(3,2,1)
plot(V1S11)
title('shelf UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(V2S11);plot(V3S11);plot(V4S11);
line([0, 40912],[0,0], 'Color', 'k');
legend('V1S','V2S','V3S','V4S')
hold off

subplot(3,2,2)
plot(V1S12)
title('slope UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(V2S12); plot(V3S12); plot(V4S12)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1S','V2S','V3S','V4S')
hold off

subplot(3,2,3)
plot(V1S21)
title('shelf LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(V2S21); plot(V3S21); plot(V4S21)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1S','V2S','V3S','V4S')
hold off

subplot(3,2,4)
plot(V1S22)
title('slope LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(V2S22); plot(V3S22); plot(V4S22)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1S','V2S','V3S','V4S')
hold off

subplot(3,2,5)
plot(V1S31)
title('shelf Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(V2S31); plot(V3S31); plot(V4S31)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1S','V2S','V3S','V4S')
hold off

subplot(3,2,6)
plot(V1S32)
title('slope Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(V2S32); plot(V3S32); plot(V4S32)
line([0, 40912],[0,0], 'Color', 'k');
legend('V1S','V2S','V3S','V4S')
hold off

%%%-----H1T,H2T,H3T,H4T------%

H1S11=H1S(1,1,:);H2S11=H2S(1,1,:);H3S11=H3S(1,1,:);H4S11=H4S(1,1,:);
H1S11=H1S11(:);H2S11=H2S11(:);H3S11=H3S11(:);H4S11=H4S11(:);

H1S21=H1S(2,1,:);H2S21=H2S(2,1,:);H3S21=H3S(2,1,:);H4S21=H4S(2,1,:);
H1S21=H1S21(:);H2S21=H2S21(:);H3S21=H3S21(:);H4S21=H4S21(:);

H1S31=H1S(3,1,:);H2S31=H2S(3,1,:);H3S31=H3S(3,1,:);H4S31=H4S(3,1,:);
H1S31=H1S31(:);H2S31=H2S31(:);H3S31=H3S31(:);H4S31=H4S31(:);

H1S12=H1S(1,2,:);H2S12=H2S(1,2,:);H3S12=H3S(1,2,:);H4S12=H4S(1,2,:);
H1S12=H1S12(:);H2S12=H2S12(:);H3S12=H3S12(:);H4S12=H4S12(:);

H1S22=H1S(2,2,:);H2S22=H2S(2,2,:);H3S22=H3S(2,2,:);H4S22=H4S(2,2,:);
H1S22=H1S22(:);H2S22=H2S22(:);H3S22=H3S22(:);H4S22=H4S22(:);

H1S32=H1S(3,2,:);H2S32=H2S(3,2,:);H3S32=H3S(3,2,:);H4S32=H4S(3,2,:);
H1S32=H1S32(:);H2S32=H2S32(:);H3S32=H3S32(:);H4S32=H4S32(:);

subplot(3,2,1)
plot(H1S11)
title('shelf UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(H2S11);plot(H3S11);plot(H4S11);
line([0, 40912],[0,0], 'Color', 'k');
legend('H1S','H2S','H3S','H4S')
hold off

subplot(3,2,2)
plot(H1S12)
title('slope UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(H2S12); plot(H3S12); plot(H4S12)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1S','H2S','H3S','H4S')
hold off

subplot(3,2,3)
plot(H1S21)
title('shelf LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(H2S21); plot(H3S21); plot(H4S21)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1S','H2S','H3S','H4S')
hold off

subplot(3,2,4)
plot(H1S22)
title('slope LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(H2S22); plot(H3S22); plot(H4S22)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1S','H2S','H3S','H4S')
hold off

subplot(3,2,5)
plot(H1S31)
title('shelf Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(H2S31); plot(H3S31); plot(H4S31)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1S','H2S','H3S','H4S')
hold off

subplot(3,2,6)
plot(H1S32)
title('slope Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(H2S32); plot(H3S32); plot(H4S32)
line([0, 40912],[0,0], 'Color', 'k');
legend('H1S','H2S','H3S','H4S')
hold off

%%%-----X1T,X2T,X3T,X4T------%

X1S11=X1S(1,1,:);X2S11=X2S(1,1,:);X3S11=X3S(1,1,:);X4S11=X4S(1,1,:);
X1S11=X1S11(:);X2S11=X2S11(:);X3S11=X3S11(:);X4S11=X4S11(:);

X1S21=X1S(2,1,:);X2S21=X2S(2,1,:);X3S21=X3S(2,1,:);X4S21=X4S(2,1,:);
X1S21=X1S21(:);X2S21=X2S21(:);X3S21=X3S21(:);X4S21=X4S21(:);

X1S31=X1S(3,1,:);X2S31=X2S(3,1,:);X3S31=X3S(3,1,:);X4S31=X4S(3,1,:);
X1S31=X1S31(:);X2S31=X2S31(:);X3S31=X3S31(:);X4S31=X4S31(:);

X1S12=X1S(1,2,:);X2S12=X2S(1,2,:);X3S12=X3S(1,2,:);X4S12=X4S(1,2,:);
X1S12=X1S12(:);X2S12=X2S12(:);X3S12=X3S12(:);X4S12=X4S12(:);

X1S22=X1S(2,2,:);X2S22=X2S(2,2,:);X3S22=X3S(2,2,:);X4S22=X4S(2,2,:);
X1S22=X1S22(:);X2S22=X2S22(:);X3S22=X3S22(:);X4S22=X4S22(:);

X1S32=X1S(3,2,:);X2S32=X2S(3,2,:);X3S32=X3S(3,2,:);X4S32=X4S(3,2,:);
X1S32=X1S32(:);X2S32=X2S32(:);X3S32=X3S32(:);X4S32=X4S32(:);

subplot(3,2,1)
plot(X1S11)
title('shelf UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(X2S11);plot(X3S11);plot(X4S11);
line([0, 40912],[0,0], 'Color', 'k');
legend('X1S','X2S','X3S','X4S')
hold off

subplot(3,2,2)
plot(X1S12)
title('slope UL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(X2S12); plot(X3S12); plot(X4S12)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1S','X2S','X3S','X4S')
hold off

subplot(3,2,3)
plot(X1S21)
title('shelf LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(X2S21); plot(X3S21); plot(X4S21)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1S','X2S','X3S','X4S')
hold off

subplot(3,2,4)
plot(X1S22)
title('slope LL')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(X2S22); plot(X3S22); plot(X4S22)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1S','X2S','X3S','X4S')
hold off

subplot(3,2,5)
plot(X1S31)
title('shelf Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(X2S31); plot(X3S31); plot(X4S31)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1S','X2S','X3S','X4S')
hold off

subplot(3,2,6)
plot(X1S32)
title('slope Dem')
xlabel('time step')
ylabel('psu.s^-1')
hold on
plot(X2S32); plot(X3S32); plot(X4S32)
line([0, 40912],[0,0], 'Color', 'k');
legend('X1S','X2S','X3S','X4S')
hold off




























