function LBbio = buildBiologicalLateralBC(LBstruct, Grd)
%BUILDBIOLOGICALLATERALBC Prepare lateral BCs for planktonic variables
%
% Each LB tracer is:
%   1) formatted using WCVIE2E_initinterpdata
%   2) interpolated to solver time step
%   3) returned as a struct of [datevec data]

nlb = numel(LBstruct);
names = cell(1, nlb);

for i = 1:nlb
    names{i} = erase(LBstruct(i).name, 'input');
end

lb = struct();

for i = 1:nlb
    name = names{i};

    A = initinterpdata( ...
        'time and boundaries', ...
        LBstruct(i).data, ...
        Grd);

    % Interpolate to solver timestep
    A.data = interp1(A.t, A.data, Grd.datatime);
    A.t    = Grd.datatime;

    lb.(name) = A;
end

% Convert time to datevec (legacy output format)
tsec = lb.(names{1}).t(:);
dnstart = datenum(Grd.start_date);
dn = dnstart + (tsec./86400);

% Assemble final output
LBbio = struct();
for i = 1:nlb
    name = names{i};
    LBbio.(name) = [dn lb.(name).data];
end

end