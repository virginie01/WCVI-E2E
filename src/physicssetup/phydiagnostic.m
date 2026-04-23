function [diagnames,flux] = phydiagnostic()

dbnames = {...
 'dtemp'     'dB/dt:temperature'       'deg C s^-1'
 'dsal'      'dB/dt:salinity'          'psu s^-1'
 };

        A.V = [...
            1     1
            2     2];
        
        A.H = [...
            1     1
            2     2];

        
        A.X = [...
            1     1
            2     2];

        
        A.CS = [...
            1     1
            2     2];

        
        A.P = [...
            1     1
            2     2];

        
        A.R = [...
            1     1
            2     2];


        A.VICC = [...
            1     1
            2     2];

        
        A.DC = [...
            1     1
            2     2];


        A.SBC = [...
            1     1
            2     2];

        
        A.CU = [...
            1     1
            2     2];

        
        A.sol_Tflx = [...
            1     1];

        
        A.srf_Tflx = [...
            1     1];
        
        flux = cell(0,3);

fld = fieldnames(A);
for ifld = 1:length(fld)
    nlnk = size(A.(fld{ifld}),1);
    
    flux = [flux; repmat(fld(ifld), nlnk, 1), num2cell(A.(fld{ifld}))];
end

nd = size(flux,1);

fluxnames = cell(nd,3);
for id = 1:nd
    fluxnames{id,1} = sprintf('%s_%02d_%02d', flux{id,:});
    fluxnames{id,2} = sprintf('%s: %02d to %02d', flux{id,:});
    fluxnames{id,3} = 'deg C/psu s^-1';
end

diagnames = [dbnames; fluxnames];


end

