function O2 = buildO2Forcing(o2_input, Grd)
%BUILDO2FORCING Prepare dissolved oxygen forcing
%
% Used when oxygen is externally prescribed (no internal O2 cycle).

O2 = initinterpdata('time and full grid', o2_input, Grd);

end