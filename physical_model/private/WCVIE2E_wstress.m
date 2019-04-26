function [taux, tauy, u10, v10]=WCVIE2E_wstress(u,v,z0)
% WSTRESS  computes wind stress using the Large and Pond (JPO, 1981) formulation.

%  USAGE: [taux,tauy,u10,v10]=wstress(u, v, z0)
%
%       u,v = nt x nx matrices. Eastward (u) and northward (v) components 
%             of wind in m/s. 
%       z0  = height of wind sensor in m (assumed to be 10 m if not
%             supplied).
%     
%       taux,tauy= nt x nx matrices. Eastward, northward wind stress 
%                  components (dynes/cm2)
%       u10,v10 =  nt x nx matrices. Equivalent wind velocity components 
%                  at 10 m.
%       
%
% Rich Signell  rsignell@usgs.gov
%
% Version 1.1 2005-04-15
% Bug fixed:  Line 49 was changed from : t=(rho*cd.*d1*1.e4)./(1+a1.*sqrt(cd));
%                                   to : t=(rho*cd.*d1*1.e4)
%             The effect of this bug was to yield slightly incorrect wind stresses for 
%             anemometer heights different from 10 m. For example, a 10 m wind at 5 m 
%             height yielded 1.806 dyn/cm2 instead of 1.699 dyn/cm2 (6.3% too high), while a 
%             10 m wind at 30 m height yielded 1.154 dyn/cm2 instead of 1.261 dyn/cm2 (8.5% too low)
%             Only the wind stresses were affected.  The output 10m winds were correct.
%             Thanks to Charlie Stock for finding this bug.

if(nargin==2)
 z0=10.;
end
%

[mu,nu] = size(u);
[mv,nv] = size(v);

if (mu~=mv) || (nu~=nv)
   error('WCVIE2E_wstress.m: u,v inputs must have the same dimensions')
end

for j = 1:size(u,2)
utmp=u(:,j);
vtmp=v(:,j);
nans=ones(length(utmp),1)*nan;
u10tmp=nans; 
v10tmp=nans;
tauxtmp=nans;
tauytmp=nans;

% identify times of zero wind.  These indices will be set to zero stress.
izeros=find((abs(utmp)+abs(vtmp))==0);
% igood=find(finite(u)& ((abs(u)+abs(v))>0));
igood = isfinite(utmp) & (abs(utmp)+abs(vtmp))>0; % KAK, finite now obsolete

utmp=utmp(igood);
vtmp=vtmp(igood);
utmp=utmp(:);
vtmp=vtmp(:);
w=utmp+1i*vtmp;
v0=abs(w);
k=.41;
rho=1.25e-3;
a1=1/k*log(z0/10.);
d1=1.e35*ones(size(utmp));
c=zeros(size(utmp));
while (max(abs(c-d1))>.01)
  c=d1;
  cd=ones(size(utmp))*1.205e-3;
  ind=find(c>11);
  if(~isempty(ind))
   cd(ind)=(0.49 + 0.065*c(ind))*1.e-3;
  end
  d1=v0./(1+a1.*sqrt(cd));
end
%t=(rho*cd.*d1*1.e4)./(1+a1.*sqrt(cd));  % wrong
t=(rho*cd.*d1*1.e4);                     % fix for version 1.1 
w10=(d1./v0).*w;
u10tmp(igood)=real(w10);
v10tmp(igood)=imag(w10);
tauxtmp(igood)=t.*u10tmp(igood); % in dynes/cm2
tauytmp(igood)=t.*v10tmp(igood); %in dynes/cm2


% set zero wind periods to zero stress (and U10,V10)
u10tmp(izeros)=zeros(size(izeros));
v10tmp(izeros)=zeros(size(izeros));
tauxtmp(izeros)=zeros(size(izeros));
tauytmp(izeros)=zeros(size(izeros));


u10(:,j)=u10tmp(:);
v10(:,j)=v10tmp(:);
taux(:,j)=tauxtmp(:);
tauy(:,j)=tauytmp(:);

end
