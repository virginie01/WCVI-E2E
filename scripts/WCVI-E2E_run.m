%% Setting up your path
% Make sure the present working directory contains the WCVI-E2E code, data,
% results, and Ecopath_Matlab subfolders. The WCVI-E2E code is spread over several 
% subfolders. To get the basics running, you need to add all subfolders to 
% the search path:

thisfile = matlab.desktop.editor.getActiveFilename;
runfolder = fileparts(thisfile);
addpath(genpath(runfolder));


%% Running physics-only simulations
%
% The simplest way to run the WCVI-E2E model is to do so without any of the 
% biological modules turned on. The code comes with some forcing datasets 
% representative of the WCVI. To run with all WCVI forcing variables, simply 
% specify an output folder name. Here, it's 'outputs_physics'. The output 
% file is the table /results/outputs_physics/archivedata.mat

physicalmodel('outputs_physics');

%% Adding biology
%
% In order to add biology, we need to turn on the biological module 'biomodel'
% via the |biofun| input, and supply any additional input parameters needed
% by the module. Output file is a table

physicalmodel('outputs_e2e','biofun',@biomodel);

% load output file
filepath = fullfile(runfolder,'results','outputs_e2e','archivedata.mat');
load(filepath);

%% define variable of interest to be plotted

% get cell array containing all variable names 
vars = file.Properties.VariableNames;

% Variable identifier to change accordingly
varname = 'PS'; 
varid  = 01;
varid = sprintf('%02d', varid);

%% Plot model outputs

% Spatial locations and titles for subplots
locations = {'ULsh', 'ULsl', 'LLsh', 'LLsl', 'DEMsh', 'DEMsl'};
titles = {'shelf UL', 'slope UL', 'shelf LL', 'slope LL', 'shelf DEM', 'slope DEM'};

% Pre-allocate variables
data = cell(1, length(locations));

% Fetch data dynamically
for i = 1:length(locations)
    var = strcat(varname, '_', locations{i});
    data{i} = file.(var);
end

% Figure 1: full simulation
figure;
for i = 1:length(data)
    subplot(3, 2, i);
    plot(data{i});
    title(titles{i});
    xlabel('Month since simulation start time');
    ylabel([varname, ' concentration (molN.m^{-3})']);
end

% Figure 2: last year only, with month labels
figure;

month_labels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

for i = 1:length(data)
    subplot(3, 2, i);

    n = length(data{i});
    last12_idx = max(1, n-11):n;
    last12_data = data{i}(last12_idx);

    x = 1:length(last12_data);
    plot(x, last12_data, '-o');

    xticks(x);
    xticklabels(month_labels(1:length(x)));

    title(titles{i});
    ylabel([varname, ' concentration (molN.m^{-3})']);
end

%% Plot derivative of variable of interest for diagnostic purposes

% Define spatial locations and titles
locations = {'ULsh', 'ULsl', 'LLsh', 'LLsl', 'DEMsh', 'DEMsl'};
titles = {'shelf UL', 'slope UL', 'shelf LL', 'slope LL', 'shelf DEM', 'slope DEM'};

% Pre-allocate cell array for data
data = cell(1, length(locations));

% Fetch data dynamically
for i = 1:length(locations)
    var = strcat(sprintf('d%s', varname), '_', locations{i});
    data{i} = file.(var);
end

% Figure 1: full simulation
figure;
for i = 1:length(data)
    subplot(3, 2, i);
    plot(data{i});
    line(xlim, [0, 0], 'Color', 'k');  % Draw horizontal line spanning full x-axis
    title(titles{i});
    xlabel('Month number');
    ylabel([varname, ' db/dt (molN.m^{-3}.s^{-1})']);
end 

% Figure 2: last year only, with month labels
figure;

month_labels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

for i = 1:length(data)
    subplot(3, 2, i);

    n = length(data{i});
    last12_idx = max(1, n-11):n;
    last12_data = data{i}(last12_idx);

    x = 1:length(last12_data);
    plot(x, last12_data, '-o');
    hold on;

    % Horizontal zero line
    line(xlim, [0, 0], 'Color', 'k');

    xticks(x);
    xticklabels(month_labels(1:length(x)));

    title([titles{i}, ' (last year)']);
    xlabel('Month');
    ylabel([varname, ' db/dt (molN.m^{-3}.s^{-1})']);

    hold off;
end

%% Plot intermediate fluxes for the corresponding variable for diagnostic purposes
%% One flux at a time

% Define spatial locations and subplot titles
locations = {'ULsh', 'ULsl', 'LLsh', 'LLsl', 'DEMsh', 'DEMsl'};
titles = {'shelf UL', 'slope UL', 'shelf LL', 'slope LL', 'shelf DEM', 'slope DEM'};

% Fluxes
fluxes = {'gpp', 'exx', 'res', 'gra', 'pre_in', 'pre_out', 'exc_in', 'exc_out', ...
    'ege_in', 'ege_out', 'mor', 'fish', 'amm', 'den', 'nit', 'vflx', 'V', ...
    'H', 'X', 'CS', 'P', 'R', 'VICC', 'DC', 'SBC', 'CU', 'sol_Tflx', 'srf_Tflx'};

month_labels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

nflx = length(fluxes);

% Iterate over each flux and produce Figure 1 = full simulation figure
for i = 1:nflx
    flx = fluxes{i};
    
    % Match flux in and flux out using regex
    flxin = ~cellfun(@isempty, regexp(vars, strcat(flx, '_\d\d_', varid, '.*')));
    flxout = ~cellfun(@isempty, regexp(vars, strcat(flx, '_', varid, '_\d\d', '.*')));
    
    % Remove overlaps (ensure `flxout` is exclusive)
    flxout_only = flxout & ~flxin;
    
    % Combine flux in and flux out
    flx_indices = find(flxin | flxout_only);
    nflx_indices = length(flx_indices);

    % Skip if no fluxes are found
    if nflx_indices == 0
        continue;
    end

    % Create a new figure for each flux
    figure;

    % Iterate over spatial locations for subplots
    for loc_idx = 1:length(locations)
        subplot(3, 2, loc_idx);
        hold on;

        loc_str = locations{loc_idx};

        % Select only fluxes relevant to this location
        flx_at_loc = flx_indices(endsWith(vars(flx_indices), loc_str));
        
        for j = 1:length(flx_at_loc)
            flxj = vars{flx_at_loc(j)};
            data = file.(flxj);
            if ismember(flxj, vars(flxout_only))
               data = -data;  % Make flux out negative
            end
            plot(data, 'DisplayName', flxj);
        end

    legend('show', 'Interpreter', 'none');
    title(titles{loc_idx});
    xlabel('Month');
    ylabel([flx, ' flux into/out of ', varname, ' (molN.m^{-3}.s^{-1})']);
    hold off;
    end
end

% Iterate over each flux and produce Figure 2 = last year figure
for i = 1:nflx
    flx = fluxes{i};
    
    % Match flux in and flux out using regex
    flxin = ~cellfun(@isempty, regexp(vars, strcat(flx, '_\d\d_', varid, '.*')));
    flxout = ~cellfun(@isempty, regexp(vars, strcat(flx, '_', varid, '_\d\d', '.*')));
    
    % Remove overlaps (ensure `flxout` is exclusive)
    flxout_only = flxout & ~flxin;
    
    % Combine flux in and flux out
    flx_indices = find(flxin | flxout_only);
    nflx_indices = length(flx_indices);

    % Skip if no fluxes are found
    if nflx_indices == 0
        continue;
    end

    % Create a new figure for each flux
    figure;

    % Iterate over spatial locations for subplots
    for loc_idx = 1:length(locations)
        subplot(3, 2, loc_idx);
        hold on;

        loc_str = locations{loc_idx};

        % Select only fluxes relevant to this location
        flx_at_loc = flx_indices(endsWith(vars(flx_indices), loc_str));
        
        for j = 1:length(flx_at_loc)
            flxj = vars{flx_at_loc(j)};
            data = file.(flxj);
            if ismember(flxj, vars(flxout_only))
               data = -data;  % Make flux out negative
            end
            n = length(data);
            last12_idx = max(1, n-11):n;
            last12_data = data(last12_idx);

            x = 1:length(last12_data);
            plot(x, last12_data, '-o', 'DisplayName', flxj);
        end

        xticks(1:length(last12_data));
        xticklabels(month_labels(1:length(last12_data)));

        legend('show', 'Interpreter', 'none');
        title(titles{loc_idx});
        xlabel('Month');
        ylabel([flx, ' flux into/out of ', varname, ' (molN.m^{-3}.s^{-1})']);
        hold off;
    end
end


%% Plot intermediate fluxes for the corresponding variable for diagnostic purposes
%% All fluxes plotted together for a given variable

% Define patterns for flux in and flux out
flux_patterns = {'...', '..._..', '....', '.', '..', '..._....'};

% Combine regex matching for fluxes in
flxin = false(size(vars)); % Initialize boolean array
for pattern = flux_patterns
    flxin = flxin | ~cellfun(@isempty, regexp(vars, strcat(pattern{1}, '_\d\d_', varid, '.*')));
end

% Combine regex matching for fluxes out
flxout_all = false(size(vars)); % Initialize boolean array
for pattern = flux_patterns
    flxout_all = flxout_all | ~cellfun(@isempty, regexp(vars, strcat(pattern{1}, '_', varid, '_\d\d', '.*')));
end

% Ensure `flxout` contains only fluxes that are exclusively flux out
flxout = flxout_all & ~flxin;

% Combine unique flux in and flux out indices
flx = flxin | flxout;
flx_indices = find(flx); % Indices of matched variables
nflx = length(flx_indices);

% Define spatial locations and subplot titles
locations = {'ULsh', 'ULsl', 'LLsh', 'LLsl', 'DEMsh', 'DEMsl'};
titles = {'shelf UL', 'slope UL', 'shelf LL', 'slope LL', 'shelf DEM', 'slope DEM'};

month_labels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

% Create figure and plot data: Full simulation figure
if nflx > 0
    figure;
    for i = 1:length(locations)
        subplot(3, 2, i);
        hold on;

        loc_str = locations{i};

        % Filter to fluxes for current location only
        flx_at_loc = flx_indices(endsWith(vars(flx_indices), loc_str));

        for j = 1:length(flx_at_loc)
            flxj = vars{flx_at_loc(j)};
            if ismember(flxj, vars(flxout)) % Check if this is flux out
                data = -file.(flxj);        % Flux out = negative
            else
                data = file.(flxj);         % Flux in = positive
            end
            h = plot(data, 'DisplayName', flxj);
            % Add flux name to data tip
            flxj_escaped = strrep(flxj, '_', '\_');

            row = dataTipTextRow('Flux', repmat({flxj_escaped}, numel(data), 1));
            h.DataTipTemplate.DataTipRows(end+1) = row;
        end

        legend('show', 'Interpreter', 'none');
        title(titles{i});
        xlabel('Month');
        ylabel(['fluxes into/out of ', varname, ' (molN.m^{-3}.s^{-1})']);
        hold off;
    end
end

% Create figure and plot data: Full simulation figure
if nflx > 0
    figure;
    for i = 1:length(locations)
        subplot(3, 2, i);
        hold on;

        loc_str = locations{i};

        % Filter to fluxes for current location only
        flx_at_loc = flx_indices(endsWith(vars(flx_indices), loc_str));

        x = []; % initialize for xticks

        for j = 1:length(flx_at_loc)
            flxj = vars{flx_at_loc(j)};
            if ismember(flxj, vars(flxout)) % Check if this is flux out
                data = -file.(flxj);        % Flux out = negative
            else
                data = file.(flxj);         % Flux in = positive
            end

            n = length(data);
            last12_idx = max(1, n-11):n;
            last12_data = data(last12_idx);

            x = 1:length(last12_data);
            h = plot(x, last12_data, '-o', 'DisplayName', flxj);

            % Add flux name to data tip
            flxj_escaped = strrep(flxj, '_', '\_');
            row = dataTipTextRow('Flux', repmat({flxj_escaped}, numel(last12_data), 1));
            h.DataTipTemplate.DataTipRows(end+1) = row;
        end

        if ~isempty(x)
            xticks(x);
            xticklabels(month_labels(1:length(x)));
        end

        legend('show', 'Interpreter', 'none');
        title(titles{i});
        xlabel('Month');
        ylabel(['fluxes into/out of ', varname, ' (molN.m^{-3}.s^{-1})']);
        hold off;
    end
end

