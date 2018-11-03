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