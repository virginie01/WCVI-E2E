function [avg, data, file] = archivedata(Grd, Arch, it, data, it0)
%ARCHIVEDATA Write physicalmodel results to file
%
% New version should be compatible with all versions of Matlab, as long as
% snctools is installed and mexnc points to the proper backend (either
% newer internal stuff or the proper mex version).
%
% Input variables:
%
%   Grd:        1 x 1 structure holding spatial and temporal grid data for
%               physical model simulation.
%
%   Arch:       1 x 1 structure holding data related to archiving.
%
%   it:         index of current model time step.
%
%   data:       nvar x 4 cell array of data to be written to file.  Column
%               1 holds the actual data values for the current model time
%               step, which must be nz x nx array.  Columns 2-4 are 
%               strings with the variables' short names, long names, and 
%               units, respectively.
%
%   it0:        index of first time step being simulated.
%
% Output variables:
%
%   avg:        nvar x 1 cell array, holding the newly-averaged values of each
%               variable.  Averages are performed over each archiving
%               period.
%
%   file:       table containing the averaged data of interest - the number
%               of rows corresponds to the number of archiving intervals, 
%               and columns hold data for each combination of data and 
%               spatial box - Variables of interest and diagnostic variables 
%               are included in the archived data


%% House-keeping constants

thisfile = mfilename('fullpath');% added
thisdir  = fileparts(thisfile); %added
repoRoot = fileparts(fileparts(thisdir)); %added

outDir   = fullfile(repoRoot, 'results', Arch.outfile); %added
matFile  = fullfile(outDir,'archivedata.mat'); % added

ndata= size(data,1);
nrow = size(Arch.dateedge(1:end-1),1); 
nz = size(Grd.z,1); nx = size(Grd.x,2);
ncol = nz.*nx.*ndata;

%% Initialize running averages the first time this function is called
if it == 1
        sz1 = cellfun(@size, data(:,1), 'uni', 0);
        avg = cellfun(@zeros, sz1, 'uni', 0);
else
        avg = Arch.avg;
end

%% Update running averages with current timestep contribution
for idata = 1:ndata
    avg{idata} = avg{idata} + Arch.fraction(it) .* data{idata};% why only one element in data{idata} if it's a matrix???Not data{idata,1}???
end

%% Create header and pre-allocated table on the very first archive call

if it == it0
    
          file = zeros(nrow,ncol);
          variablenames = cell(1,ncol);
          for i=1:ndata
              variablenames(6*(i-1)+1:6*i) = {strcat(data{i,2},'_ULsh'),strcat(data{i,2},'_LLsh'),strcat(data{i,2},'_DEMsh'),strcat(data{i,2},'_ULsl'),strcat(data{i,2},'_LLsl'),strcat(data{i,2},'_DEMsl')};
          end
          
          file = array2table(file,'VariableNames',variablenames);

          if ~exist(outDir,'dir'), mkdir(outDir); end % added
          save(matFile, 'file'); % added
        
else
    
     file = Arch.file;

end

%% End-of-window: plug averages & append to disk

if Arch.islast(it)
    
    vec = zeros(1,ncol);
    for iv = 1:length(avg)
        vectmp = reshape(avg{iv},1,[]);
        id1 = (6*(iv-1))+1;
        id2 = (6*iv);
        vec(id1:id2) = vectmp;
    end
    
     file{Arch.bin(it),:} = vec; 

     % --- append one row to .mat file without rewriting the header -------------
     save(matFile,'file','-append'); % added
    
     % reset running means
    
     avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0); 
    
end



    



