%% Define the paths for the project

% type in the root path for the data
root_path = '';

% list of paths
path_list = {'clusters','param','rois',...
    'registration','classifier','reference','external'};
% allocate the structure
paths = struct([]);

% add the root
paths(1).main_path = root_path;
% get the number of paths
path_number = length(path_list);

% for all the paths on the list
for path_var = 1:path_number
    % assemble the field name
    field = strcat(path_list{path_var},'_path');
    % insert the path in the structure
    paths(1).(field) = ...
        fullfile(root_path,'Analysis',path_list{path_var},'\');
end

% add the AF10 and Tectum colors
paths(1).afOT_colors = [0.196078431372549,0.345098039215686,0.400000000000000;...
    0.713174463000000,0.213888935296351,0.476258604429344];