function LBbio = diapauseLBbio(LBbio, B, names, A)
% LBbiodiapause enables to calculate the concentration of ZL1 and ZL2 at the lateral boundaries when
% diapause is on. ZL1 and ZL2 are determined based on the total concentration of ZL in a given 
% source and the proportion of ZL1 and ZL2 in the destination box

% INPUTS
%
% it:      current time step
% 
% Grd:     Structure holding spatial and temporal grid data for WCVIE2E_physicalmodel simulations
%
% LBbio:   Nested structure. First level = "zooplanktonic" group (i.e. mixed). Second level:
%          
%          t:    Grd.time(it):In.datadt:Grd.time(it+1), time (seconds from sim start time)
%          o:    1 x 6 array specifying the boundary
%          data: length(t) x length(o) array holding the data for the current time steps
%
% B:       nz x nx x nbsv array. current concentrations of biological variables.
%
% A:       structure containing biological parameters.
%
%
% OUTPUTS
%
% LBbiodiapause: new updated LBbio with appropriate calculated values for
%                groups ZL, ZL1 and ZL2 from time Grd.time(it) to Grd.time(it+1)


ZL = names{A.idx.zl,2};
ZL1 = names{A.idx.zl1,2};
ZL2 = names{A.idx.zl2,2};

n = length(LBbio.(ZL).t);

for i = 1:n
    
    % source 1 = open ocean upper layer. Destination box = slope upper layer

    LBbio.(ZL1).data(i,1) = LBbio.(ZL).data(i,1).* ...
        (B(1,2,A.idx.zl1)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL2).data(i,1) = LBbio.(ZL).data(i,1).* ...
        (B(1,2,A.idx.zl2)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL).data(i,1) = 0;
    
    % source 2 = open ocean lower layer. Destination box = slope lower layer
    
    LBbio.(ZL1).data(i,2) = LBbio.(ZL).data(i,2).* ...
        (B(2,2,A.idx.zl1)./(B(2,2,A.idx.zl1)+B(2,2,A.idx.zl2)));
    LBbio.(ZL2).data(i,2) = LBbio.(ZL).data(i,2).* ...
        (B(2,2,A.idx.zl2)./(B(2,2,A.idx.zl1)+B(2,2,A.idx.zl2)));
    LBbio.(ZL).data(i,2) = 0;
    
    % source 3 = rain over the shelf. Destination box = shelf upper layer
    
    LBbio.(ZL1).data(i,3) = LBbio.(ZL).data(i,3).* ...
        (B(1,1,A.idx.zl1)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,3) = LBbio.(ZL).data(i,3).* ...
        (B(1,1,A.idx.zl2)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL).data(i,3) = 0;
 
    % source 4 = rain over the slope. Destination box = slope upper layer

    LBbio.(ZL1).data(i,4) = LBbio.(ZL).data(i,4).* ...
        (B(1,2,A.idx.zl1)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL2).data(i,4) = LBbio.(ZL).data(i,4).* ...
        (B(1,2,A.idx.zl2)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL).data(i,4) = 0;
    
    % source 5 = freshwater from run-offs. Destination box = shelf upper layer
    
    LBbio.(ZL1).data(i,5) = LBbio.(ZL).data(i,5).* ...
        (B(1,1,A.idx.zl1)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,5) = LBbio.(ZL).data(i,5).* ...
        (B(1,1,A.idx.zl2)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL).data(i,5) = 0;
   
    % source 6 = VICC. Destination box = shelf upper layer

    LBbio.(ZL1).data(i,6) = LBbio.(ZL).data(i,6).* ...
        (B(1,1,A.idx.zl1)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,6) = LBbio.(ZL).data(i,6).* ...
        (B(1,1,A.idx.zl2)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL).data(i,6) = 0;
    
    % source 7 = sbc shelf ul. Destination box = shelf UL

    LBbio.(ZL1).data(i,7) = LBbio.(ZL).data(i,7).* ...
        (B(1,1,A.idx.zl1)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,7) = LBbio.(ZL).data(i,7).* ...
        (B(1,1,A.idx.zl2)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL).data(i,7) = 0;
    
    % source 8 = sbc shelf ll. Destination box = shelf lower layer

    LBbio.(ZL1).data(i,8) = LBbio.(ZL).data(i,8).* ...
        (B(2,1,A.idx.zl1)./(B(2,1,A.idx.zl1)+B(2,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,8) = LBbio.(ZL).data(i,8).* ...
        (B(2,1,A.idx.zl2)./(B(2,1,A.idx.zl1)+B(2,1,A.idx.zl2)));
    LBbio.(ZL).data(i,8) = 0;
    
    % source 9 = sbc slope ul. Destination box = slope upper layer

    LBbio.(ZL1).data(i,9) = LBbio.(ZL).data(i,9).* ...
        (B(1,2,A.idx.zl1)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL2).data(i,9) = LBbio.(ZL).data(i,9).* ...
        (B(1,2,A.idx.zl2)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL).data(i,9) = 0;
    
    % source 10 = dc shelf ul. Destination box = shelf upper layer

    LBbio.(ZL1).data(i,10) = LBbio.(ZL).data(i,10).* ...
        (B(1,1,A.idx.zl1)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,10) = LBbio.(ZL).data(i,10).* ...
        (B(1,1,A.idx.zl2)./(B(1,1,A.idx.zl1)+B(1,1,A.idx.zl2)));
    LBbio.(ZL).data(i,10) = 0;
    
    % source 11 = dc shelf ll. Destination box = shelf lower layer

    LBbio.(ZL1).data(i,11) = LBbio.(ZL).data(i,11).* ...
        (B(2,1,A.idx.zl1)./(B(2,1,A.idx.zl1)+B(2,1,A.idx.zl2)));
    LBbio.(ZL2).data(i,11) = LBbio.(ZL).data(i,11).* ...
        (B(2,1,A.idx.zl2)./(B(2,1,A.idx.zl1)+B(2,1,A.idx.zl2)));
    LBbio.(ZL).data(i,11) = 0;
    
    % source 12 = dc slope ul. Destination box = slope upper layer

    LBbio.(ZL1).data(i,12) = LBbio.(ZL).data(i,12).* ...
        (B(1,2,A.idx.zl1)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL2).data(i,12) = LBbio.(ZL).data(i,12).* ...
        (B(1,2,A.idx.zl2)./(B(1,2,A.idx.zl1)+B(1,2,A.idx.zl2)));
    LBbio.(ZL).data(i,12) = 0;
    
    % source 13 = cu slope ll. Destination box = slope lower layer

    LBbio.(ZL1).data(i,13) = LBbio.(ZL).data(i,13).* ...
        (B(2,2,A.idx.zl1)./(B(2,2,A.idx.zl1)+B(2,2,A.idx.zl2)));
    LBbio.(ZL2).data(i,13) = LBbio.(ZL).data(i,13).* ...
        (B(2,2,A.idx.zl2)./(B(2,2,A.idx.zl1)+B(2,2,A.idx.zl2)));
    LBbio.(ZL).data(i,13) = 0;








    
end

end

