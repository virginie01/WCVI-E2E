function Grd = buildGrid(In)
%BUILDGRID Construct spatial and temporal grids for WCVI-E2E model

% --- Vertical grid ---
In.dz = In.dz(:,:);
Grd.zp = [zeros(1,size(In.dz,2)); cumsum(In.dz,1)];
Grd.z  = (Grd.zp(1:end-1,:) + Grd.zp(2:end,:)) ./ 2;
Grd.nz = size(Grd.z,1);
Grd.boxz = {'Upper Layer'; 'Lower Layer'; 'Demersal'};

% --- Horizontal grid ---
In.dx = In.dx(:,:);
Grd.xp = cumsum(In.dx,2);
Grd.xp = [zeros(size(Grd.xp,1),1), Grd.xp];
Grd.x  = (Grd.xp(:,1:end-1) + Grd.xp(:,2:end)) ./ 2;
Grd.nx = size(Grd.x,2);
Grd.boxx = {'Shelf','Slope'};

% --- Start/end dates ---
if isequal(size(In.syear),[1 6])
    Grd.start_date = In.syear;
else
    Grd.start_date = [In.syear 1 1 0 0 0];
end

if isequal(size(In.eyear),[1 6])
    Grd.end_date = In.eyear;
else
    Grd.end_date = [In.eyear 12 31 24 0 0];
end

dnstart = datenum(Grd.start_date);
dnend   = datenum(Grd.end_date);

% --- Time grids ---
Grd.tmax = (dnend - dnstart) * 86400;
Grd.nt   = floor(Grd.tmax / In.dt);
Grd.time = (0:Grd.nt).' * In.dt;

Grd.datant   = floor(Grd.tmax / In.datadt);
Grd.datatime = (0:Grd.datant).' .* In.datadt;
end