%% 0) Registration analysis
clearvars
close all

% define the figure path
Paths
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).clusters_path;
registration_path = paths(1).registration_path;

data = load_clusters(cluster_path);

% define the stimulus labels
% get the dataset name
stim_name = data.name;
%scan for the p17b
if contains(stim_name, 'p17b')
    %if it's p17b
%     stim_labels = {'Red','Green','Blue','UV'};
    stim_labels = {'R','G','B','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'R CK','UV CK','R GR','UV GR','R FL','UV FL'};
end
% get the stimulus number
stim_num = data.stim_num;
% Load transformation matrix

%define the main loading path
tmat_path = fullfile(registration_path,'Registration_info',data.name);

% get the list of files
file_list = dir(tmat_path);
file_list = file_list(3:end);

% get the number of files
num_files = length(file_list);

%allocate memory for the affine objects
aff_cell = cell(num_files,1);

%go file by file extracting the matrix
%allocate memory for the matrices
x_mat = zeros(5,3,num_files);
%for all the files
for files = 1:num_files
    %get the file name
    target_folder = file_list(files).name;

    %load the actual matrix
    temp_cell = importdata(fullfile(tmat_path,target_folder,'registration'));
    %and parse it
    %allocate memory for the parsed data
    parse_mat = zeros(5,3);
    %for all the relevant lines
    for plines = 6:10
        %split the string via spaces
        temp_str = strsplit(temp_cell{plines},' ');
        %convert the last 3 to numbers and store
        parse_mat(plines-5,:) = [str2double(temp_str{2});str2double(temp_str{3});str2double(temp_str{4})];
    end
    %store the parse matrix in the main matrix
    x_mat(:,:,files) = parse_mat;

end   
% Turn the transformation matrix into an affine matrix

%allocate memory for the transformation matrices
trans_mats = cell(num_files,1);

%for all the files
for files = 1:num_files
    %turn the cmtk parameter matrix into an affine matrix and store
    trans_mats{files} = cmtkparams2affmat(x_mat(:,:,files));

end

%get the corresponding affine object
%for all the files
for files = 1:num_files
    aff_cell{files} = affine3d(trans_mats{files});
end
% Load the reference brain

% define the ref brain depending on the file
if contains(data(1).name,{'Syn','syn'})
    ref_name = 'refisl2cut2.nrrd.tif';
    full_ref_name = 'Ref20131120pt14pl2.nrrd.tif';
else
    ref_name = 'refcutblursub.nrrd.tif';
    full_ref_name = 'Ref20131120pt14pl2.nrrd.tif';

end
% assemble the reference path
ref_path = fullfile(registration_path,'Reference_brains',ref_name);
full_ref_path = fullfile(registration_path,'Reference_brains',full_ref_name);
%allocate memory for the stack
ref_info = imfinfo(ref_path);
ref_stack = zeros(ref_info(1).Height,ref_info(1).Width,size(ref_info,1));

full_ref_info = imfinfo(full_ref_path);
full_ref_stack = zeros(full_ref_info(1).Height,full_ref_info(1).Width,size(full_ref_info,1));
%load the stack
%for all the z sections
for z = 1:size(ref_info,1)
    ref_stack(:,:,z) = imread(ref_path,z);
end

%for all the z sections
for z = 1:size(full_ref_info,1)
    full_ref_stack(:,:,z) = imread(full_ref_path,z);
end

%get the dimensions of the ref brain to allocate for the registered brains
full_ref_dim = size(full_ref_stack);
% Load the pre-reg tif (from pre-reg NRRD file)

%define the path
prereg_search = fullfile(registration_path,'Pre_registration_brains',data.name,'*.tif');
%get the tif files in the directory
tif_list = dir(prereg_search);
%get the number of tifs
tif_num = length(tif_list);
%allocate memory to store the prereg coords
prereg_coord = cell(tif_num,2);
% also for the im_info
im_info_cell = cell(tif_num,1);
% get the path
prereg_path = tif_list(1).folder;
%for all the files
for tifs = 1:tif_num
    
    %define the filename
    prereg_name = tif_list(tifs).name;
    
    %allocate memory for the stack
    im_info = imfinfo(fullfile(prereg_path,prereg_name));
    prereg_stack = zeros(im_info(1).Height,im_info(1).Width,size(im_info,1));
    %load the stack
    %for all the z sections
    for z = 1:size(im_info,1)
        prereg_stack(:,:,z) = imread(fullfile(prereg_path,prereg_name),z);
    end

    %also extract its coordinates
    ind_coord = find(prereg_stack);
    [x,y,z] = ind2sub(size(prereg_stack),ind_coord);
    prereg_coord{tifs,1} = prereg_stack;
    prereg_coord{tifs,2} = cat(2,x,y,z);
    % store the stack information
    im_info_cell{tifs} = im_info;
end
% Load the labels

% define the path to the file
labels_file = fullfile(registration_path,'Labels_info','MaskDatabase.mat');
% load the labels file
labels_data = load(labels_file);

% define which dataset to load depending on the reference
if contains(data.name,{'Syn','syn'})
    field_list = {'AF4','AF5','AF6','AF7','AF8','AF9','Tecum Neuropil','Periventriculare'};
    coordinate_conversion = [274,0,0];
else
    field_list = {'Tectum Stratum','Tecum Neuropil','Pretectum','Habenula','Cerebellum'};
    coordinate_conversion = [272,74,49];
end

% get the number of fields
field_number = length(field_list);
% allocate memory for the fields
label_cell = cell(field_number,2);
label_stack = zeros(labels_data.height,labels_data.width,labels_data.Zs);
label_outline = zeros(labels_data.height,labels_data.width,labels_data.Zs);
% get the field numbers for these names
for field = 1:length(field_list)
    % get the index with the first name match (so as to not get
    % subdivisions)
    idx_vector = find(contains(labels_data.MaskDatabaseNames,field_list{field}),1);
    % store the map and the name in the cell
    label_stack = label_stack + reshape(full(labels_data.MaskDatabase(:,idx_vector)),...
        labels_data.height,labels_data.width,labels_data.Zs).*field;
    % also store the outlines
    label_outline = label_outline + reshape(full(labels_data.MaskDatabaseOutlines(:,idx_vector)),...
        labels_data.height,labels_data.width,labels_data.Zs).*field;
    
    label_cell{field,2} = labels_data.MaskDatabaseNames{idx_vector};
end
% Correct for fish of origin (p17b_syngc6s in particular)

% get the fish of origin
fish_ori = data.fish_ori;
% correct the data in case there were multiple experiments per fish
% get the number of experiments
num_experiment = sum(diff(fish_ori(:,2))<0)+1;
num_fish = length(unique(fish_ori(:,1)));
% if the numbers are different, correct
if num_fish < num_experiment
    % allocate memory for the new fish_ori
    fish_ori_corrected = zeros(size(fish_ori));
    fish_ori_corrected(:,2) = fish_ori(:,2);
    % find the junctions
    junctions = [0;find(diff(fish_ori(:,2))<0);size(fish_ori,1)];

    % for all the junctions
    for junc = 1:length(junctions)-1
        fish_ori_corrected(junctions(junc)+1:junctions(junc+1),1) = junc;
    end
    % replace the old fish_ori with the new one for registration only
    fish_ori = fish_ori_corrected;
end
% Reformat the data

fishave_cell = reformat_fish(data,fish_ori,num_files,im_info_cell,aff_cell,ref_info,'seeds');
% Determine the region for each seed from the registration

% allocate memory for the registred anatomy
registered_anatomy = cell(num_files,1);
label_copy = zeros(size(label_stack));
% for all the experiments
for files = 1:num_files
    % load the coordinates
    current_coord = fishave_cell{files,1};
    % correct the coordinates for the different references
    current_coord = round(current_coord + coordinate_conversion);
    
    % record the region for each seed
    % allocate memory for the seeds
    seed_region = zeros(size(current_coord,1),1);
    % for all the seeds
    for seeds = 1:size(current_coord,1)  

        if current_coord(seeds,1) > size(label_stack,1) || ...
                current_coord(seeds,2) > size(label_stack,2) || ...
                current_coord(seeds,3) > size(label_stack,3)
            continue
        else

            seed_region(seeds) = label_stack(current_coord(seeds,1),...
                current_coord(seeds,2),current_coord(seeds,3));
        end
    end
    % store the result
    registered_anatomy{files} = seed_region;
end

% turn the anatomy into a vector
registered_anatomy = cat(1,registered_anatomy{:});
% Get the coordinate vector for each map

% combine the coordinate vectors
coord = round(cat(1,fishave_cell{:,1}))+coordinate_conversion;
% take the negative seeds out (also in the index)

% get a vector to keep the seeds only within the volume
selection_vector = coord(:,1)<full_ref_dim(1) & coord(:,2)<full_ref_dim(2) &...
       coord(:,3)<full_ref_dim(3);
selection_vector = selection_vector & sum(coord<1,2)==0;

% get the seed indexes corresponding to each pixel
indexes = cat(1,fishave_cell{:,2});
indexes = indexes(selection_vector);
coord = coord(selection_vector,:);
registered_anatomy = registered_anatomy(selection_vector);

% turn the coordinates into linear indexes
coord = sub2ind(full_ref_dim,coord(:,1),coord(:,2),coord(:,3));
%% 1) Plot maps of the max response
    
% get the reshaped activity
conc_trace = reshape(data.conc_trace,[],data.time_num,data.stim_num);
% take only the stimulation time
conc_trace = conc_trace(:,21:60,:);
% take the absolute average
conc_trace = squeeze(mean(abs(conc_trace),2));
% allocate memory for the percentile threshold
perc_threshold = ones(size(conc_trace,1),stim_num)==1;

% get the max gain for each seed
[max_gain,max_idx] = max(conc_trace,[],2);
% get the number of stimuli
stim_num = data.stim_num;
% allocate memory to store each map
gain_maps = cell(stim_num,1);
for color = 1:stim_num
    gain_maps{color} = zeros(full_ref_dim);
    % get the percentile threshold
    perc_threshold(:,color) = conc_trace(:,color) > prctile(conc_trace(:,color),75);
    
end
% run through all the seeds
for seeds = 1:size(conc_trace,1)
    
    % get the indexes for the corresponding seeds and anatomy
    index_vector = coord(indexes==seeds&registered_anatomy>0);

    % for each color
    for color = 1:stim_num
        % include the seed only if it passes the threshold for that color
        if perc_threshold(seeds,color) == 0
            continue
        end
        % assemble the gain map based on the adding the gains for each
        % voxel across animals
        gain_maps{color}(index_vector) = ...
            gain_maps{color}(index_vector) + abs(conc_trace(seeds,color));
        
    end
    
end
% allocate memory to save the projections
color_cell = cell(stim_num,1);


% trim the label stack
if contains(data.name,'p17b')
    half_label = max(label_outline(140:760,:,:)>0,[],3);
else
    half_label = max(label_outline(300:700,110:511,:)>0,[],3);
end

% generate the maps
for color = 1:stim_num
    % normalize the full ref stack as background
    curr_stack = sum(gain_maps{color},3);
    curr_blank = zeros(size(curr_stack));
    
    if contains(data.name,'p17b')
        switch color
            case 1
                target_stack = cat(4,curr_stack,curr_blank,curr_blank);
            case 2
                target_stack = cat(4,curr_blank,curr_stack,curr_blank);
            case 3
                target_stack = cat(4,curr_blank,curr_blank,curr_stack);
            case 4
                target_stack = cat(4,curr_stack,curr_blank,curr_stack);
        end
    else
        switch color
            case {1,3,5}
                target_stack = cat(4,curr_stack,curr_blank,curr_blank);
            case {2,4,6}
                target_stack = cat(4,curr_stack,curr_blank,curr_stack);
        end
    end
    
    % take the anterior half
    if contains(data.name,'p17b')
        half_stack = target_stack(140:760,:,:,:);
    else
        half_stack = target_stack(300:700,110:511,:);
    end
    
    % produce a max intensity projection
    max_projection = imgaussfilt(squeeze(half_stack),5);
    color_cell{color} = max_projection;
    
    
end
clear gain_maps


% for all the clusters
for color = 1:stim_num
    figure
   
    % create the image
    I = squeeze(color_cell{color});
    imagesc(normr_1(I,1)+half_label.*0.2)
    set(gca,'XTick',[],'YTick',[])
    axis equal
    box off

end
%% 2) Calculate the clusters one by one
close all
% get the clusters
idx_clu = data.idx_clu;
% get the number of clusters
clu_num = data.clu_num;

% get a color map for the clusters
cluster_color = distinguishable_colors(clu_num,[0 0 0;1 1 1]);
% allocate en empty stack
norm_stack = zeros([full_ref_dim,3]);

% trim the label stack
if contains(data.name,'p17b')
    half_label = max(label_outline(140:760,:,:)>0,[],3);
else
    half_label = max(label_outline(300:700,:,:)>0,[],3);
end

% allocate memory for the max projections
cluster_cell = cell(clu_num,1);
% for all the clusters
for clu = 1:clu_num
    fprintf(strjoin({'Current cluster',num2str(clu),'\r\n'},'_'))
    % allocate memory for the stack
    temp_stack = norm_stack;

    % get the seeds for this cluster
    seed_list = find(idx_clu==clu);

    % get the color of the cluster
    seed_color = cluster_color(clu,:);
    
    index_vector = coord(ismember(indexes,seed_list)&(registered_anatomy>0));
    
    [x,y,z] = ind2sub(full_ref_dim,index_vector);
    index_vector_r = sub2ind([full_ref_dim,3],x,y,z,ones(size(x)));
    index_vector_g = sub2ind([full_ref_dim,3],x,y,z,ones(size(x)).*2);
    index_vector_b = sub2ind([full_ref_dim,3],x,y,z,ones(size(x)).*3);
    
    temp_stack(index_vector_r) = seed_color(1);
    temp_stack(index_vector_g) = seed_color(2);
    temp_stack(index_vector_b) = seed_color(3);

    % take the anterior half (protocol dependent)
    if contains(data.name,'p17b')
        half_stack = temp_stack(140:760,:,:,:);
    else
        half_stack = temp_stack(300:700,:,:,:);
    end

    % AND the channels
    max_projection = any(half_stack,4);
    max_projection = sum(max_projection,3);
    % store the projection
    cluster_cell{clu} = max_projection;
end
% Plot the cluster projections one by one 
close all

% for all the clusters
for clu = 1:clu_num
    figure
    % create the image
    I = squeeze(cluster_cell{clu});
    % gaussian blur
    I = normr_1(imgaussfilt(I, 20),1)+half_label.*0.5;
    
    imagesc(I)
    colormap(magma)
    set(gca,'XTick',[],'YTick',[])
    axis square
    box off
end
%% 3) Calculate single fish maps for the AP and DV projections

if contains(data.name,'p17b')

% get the single reps
all_single_reps = data.single_reps;
% allocate memory for the matrices
single_matrices = zeros([full_ref_dim,stim_num,num_fish]);
% get the percentage threshold for all the seeds

% for all the fish
for fish = 1:num_fish
    fprintf(strjoin({'Current fish:',num2str(fish),'\r\n'},' '))
    % get the traces for this fish
    single_reps = all_single_reps(data.fish_ori(:,1)==fish,:);
    % get a seed map
    seed_map = find(data.fish_ori(:,1)==fish)';
    % determine the filtering threshold
    all_reps = abs(reshape(single_reps,size(single_reps,1),data.time_num,stim_num,[]));
    % get the number of reps
    rep_num = size(all_reps,4);

    rep_maps = zeros([full_ref_dim,stim_num,rep_num]);
    % fir all the reps
    for reps = 1:rep_num
        % get the current rep
        conc_trace = all_reps(:,:,:,reps);

        % take only the stimulation time
        conc_trace = conc_trace(:,21:60,:);
        % take the absolute average
        conc_trace = squeeze(mean(abs(conc_trace),2));
        % calculate the maps
        maps_cell = property_map(coord,indexes,registered_anatomy,...
            conc_trace,stim_num,full_ref_dim,seed_map,perc_threshold);

        % save in the rep map concatenated by stimulus
        rep_maps(:,:,:,:,reps) = cat(4,maps_cell{:});
    end

    % sum the rep maps into a single matrix produce the sectional maps
    single_matrices(:,:,:,:,fish) = sum(rep_maps~=0,5);
    
end

clear all_reps

% Calculate the projections

    % allocate memory to save the AP and DV projections
    apdv_cell = cell(stim_num,3,2);
    % allocate memory for the individual fish for stats
    fish_cell = cell(stim_num,3);
    
    % for all stimuli
    for stim = 1:stim_num
        
        % take the anterior half (protocol dependent)
        if contains(data.name,'p17b')
            half_stack = squeeze(single_matrices(140:760,1:310,:,stim,:));
            colors = [1 0 0;0 1 0;0 0 1;1 0 1];
        else
            half_stack = squeeze(single_matrices(300:700,1:310,:,stim,:));
            colors = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
        end
        for dims = 1:3
            switch dims
                case 1
                    temp_stack = half_stack;
                    rotation = 0;
                case 2
                    temp_stack = permute(half_stack,[2 1 3 4]);
                    rotation = 0;
                case 3
                    temp_stack = permute(half_stack,[3 1 2 4]);
                    rotation = 0;
            end
            % calculate the profile
            % collapse in the third dimension
            fish_profiles = squeeze(sum(temp_stack,3));
            % rotate to match the retinotopic axis
            % allocate memory for the rotations
            rotation_cell = cell(num_fish,1);
            % for all the fish
            for fish = 1:num_fish
                rotation_cell{fish} = imrotate(fish_profiles(:,:,fish),rotation);
            end
            % collapse the last dimension
            fish_profiles = squeeze(sum(cat(3,rotation_cell{:}),2));
            
            % get the average and std across fish
            apdv_cell{stim,dims,1} = mean(fish_profiles,2);
            apdv_cell{stim,dims,2} = std(fish_profiles,0,2)./sqrt(num_fish);
            % store the individual fish for stats
            fish_cell{stim,dims} = fish_profiles;
        end
        
    end

% Compare the profiles statistically
    % define the number of bins to use
    nbins = 5;
    % allocate memory for the splits
    split_cell = cell(3,2);
    for dims = 1:3
        figure
        % allocate memory for the comparisons
        comparison_matrix = zeros(stim_num,stim_num,nbins);
        switch dims
            case 1
                x_label = {'A','P'};
                start = 171;
                stop = 531;
            case 2
                x_label = {'L','M'};
                start = 112;
                stop = 310;
            case 3
                x_label = {'D','V'};
                start = 67;
                stop = 124;
        end
        
        % define the splits
        delta = stop-start+1;
        split_int = round([1:delta/nbins:delta;...
            delta/nbins:delta/nbins:delta+1]);
        % save for later
        split_cell{dims,1} = split_int;
        % for all the bins
        for bins = 1:nbins
            % assemble the matrix for plotting
            cat_matrix = cat(3,fish_cell{:,dims});
            % trim to the current bin
            cat_matrix = cat_matrix(start:stop,:,:);
            cat_matrix = cat_matrix(split_int(1,bins):split_int(2,bins),:,:);
            
            % normalize the matrices
            matrix_max = max(cat_matrix(:));
            cat_matrix = cat_matrix./matrix_max;
            
            % get the list of combinations
            list_comb = nchoosek(1:stim_num,2);
            % get the number of combinations
            number_comb = size(list_comb,1);
            
            % for all the combinations
            for combo = 1:number_comb
                
                % get the corresponding profiles
                prof_1 = reshape(cat_matrix(:,:,list_comb(combo,1)),[],1);
                prof_2 = reshape(cat_matrix(:,:,list_comb(combo,2)),[],1);
                
                % run the comparison
                comparison_matrix(list_comb(combo,1),list_comb(combo,2),bins) = signrank(prof_1,prof_2).*number_comb;
                comparison_matrix(list_comb(combo,2),list_comb(combo,1),bins) = signrank(prof_1,prof_2).*number_comb;
            end
        end
        
        comparison_matrix = comparison_matrix<0.05;
        % store the matrix
        split_cell{dims,2} = comparison_matrix;
        % plot the matrix
        for bins  = 1:nbins
            subplot(round(sqrt(nbins)),ceil(sqrt(nbins)),bins)
            imagesc(comparison_matrix(:,:,bins))
            colormap(jet)
        end
    end
    % Plot the gradients
    
    close all
    
    % define the projection labels
    projection_labels = {'AP','LM','DV'};
    
    
    for dims = 1:3
        
        switch dims
            case 1
                x_label = {'A','P'};
                start = 171;
                stop = 531;
            case 2
                x_label = {'L','M'};
                start = 112;
                stop = 310;
            case 3
                x_label = {'D','V'};
                start = 67;
                stop = 124;
        end
        % assemble the matrix for plotting
        ave_matrix = cat(2,apdv_cell{:,dims,1});
        sem_matrix = cat(2,apdv_cell{:,dims,2});
        
        ave_matrix = ave_matrix(start:stop,:);
        sem_matrix = sem_matrix(start:stop,:);
        
        % normalize the matrices
        matrix_max = max(ave_matrix(:));
        ave_matrix = ave_matrix./matrix_max;
        sem_matrix = sem_matrix./matrix_max;
        
        % define the offset
        offset = 2.2;
        % define the offset for the confidence lines
        mini_offset = -0.2;
        
        figure
        % for all the stimuli
        for stim = 1:stim_num
            % get the effective offset
            eff_offset = stim_num*offset-offset*(stim-1);
            shadedErrorBar(1:length(ave_matrix),ave_matrix(:,stim)+eff_offset,sem_matrix(:,stim),...
                {'Color',colors(stim,:)},0)
            hold on

            % get the splits
            split_int = split_cell{dims,1};
            % reset the offset
            mini_count = mini_offset;
            % for the remaining stim
            for stim2 = 1:stim_num
                % if it's the current stim, skip
                if stim == stim2
                    continue
                end
                % for all the bins
                for bins = 1:nbins
                    p_vector = split_cell{dims,2}(stim,stim2,bins);
                    
                    % define the color
                    if p_vector
                        color = colors(stim2,:);
                    else
                        continue
                    end
                    % get the line coordinates
                    x = split_int(1,bins):split_int(2,bins);
                    y = eff_offset.*ones(size(x))+mini_count;
                    
                    plot(x,y,'Color',color)
                    
                end
                % update the offset
                mini_count = mini_count + mini_offset;
            end
            
        end
        
        set(gca,'YTick',[],'Visible','off')
        
    end
end
%% 4) Calculate the overlap between stimuli and compare to reps

% get the single reps
all_single_reps = data.single_reps;

% allocate memory for the matrices
single_matrices = cell(num_fish,2);
% for all the fish
for fish = 1:num_fish
    fprintf(strjoin({'Current fish:',num2str(fish),'\r\n'},' '))
    % get the traces for this fish
    single_reps = all_single_reps(data.fish_ori(:,1)==fish,:);
    % get a seed map
    seed_map = find(data.fish_ori(:,1)==fish)';
    % determine the filtering threshold
    all_reps = abs(reshape(single_reps,size(single_reps,1),data.time_num,stim_num,[]));
    % get the number of reps
    rep_num = size(all_reps,4);
    % allocate memory to store the maps per rep
    rep_maps = cell(rep_num,1);
    % fir all the reps
    for reps = 1:rep_num
        % get the current rep
        conc_trace = all_reps(:,:,:,reps);

        % take only the stimulation time
        conc_trace = conc_trace(:,21:60,:);
        % take the absolute average
        conc_trace = squeeze(mean(abs(conc_trace),2));
        % calculate the maps
        maps_cell = property_map(coord,indexes,registered_anatomy,conc_trace,stim_num,full_ref_dim,seed_map);

        % take the absolute average
        ROI_threshold = prctile(conc_trace,75,1);

        % threshold them

        % for all the maps
        for stim = 1:stim_num
            % apply the threshold and save
            maps_cell{stim} = abs(maps_cell{stim});
            maps_cell{stim} = maps_cell{stim}(:)>ROI_threshold(stim);
        end
        % save in the rep map concatenated by stimulus
        rep_maps{reps} = cat(2,maps_cell{:});
    end

    % concatenate the rep maps into a single matrix
    rep_maps = cat(3,rep_maps{:});
    %% Calculate the overlap matrices

    % determine the number of stimulus pairs
    stim_combinations = nchoosek(1:stim_num,2);
    number_combinations = size(stim_combinations,1);

    % allocate memory for the overlaps
    overlap_matrix = zeros(stim_num,stim_num,rep_num);
    % for all the reps
    for reps = 1:rep_num
        % for all the combinations
        for combs = 1:number_combinations
            % get the maps
            map1 = rep_maps(:,stim_combinations(combs,1),reps);
            map2 = rep_maps(:,stim_combinations(combs,2),reps);

            % quantify the overlap
            overlap_matrix(stim_combinations(combs,1),stim_combinations(combs,2),reps) = sum(sum([map1,map2],2)>1);

        end
    end
    
    % save the matrix
    single_matrices{fish,1} = overlap_matrix;
    %% Calculate the control matrix with the reps for each stimulus

    % allocate memory to store the overlaps between reps
    control_overlap = zeros(stim_num,1);
    % get the combinations for the reps
    control_combs = nchoosek(1:rep_num,2);
    control_num = size(control_combs,1);
    % for all the stimuli
    for stim = 1:stim_num
        % get the maps for this stim
        stim_maps = squeeze(rep_maps(:,stim,:));
        % for all the combinations
        for combs = 1:control_num
            % select the corresponding maps
            map1 = stim_maps(:,control_combs(combs,1));
            map2 = stim_maps(:,control_combs(combs,2));
            % calculate the overlap
            control_overlap(stim) = control_overlap(stim) + sum(sum([map1,map2],2)>1)./control_num;
        end
    end

    % save the control matrix
    single_matrices{fish,2} = control_overlap;
end

clear all_reps
% Calculate significance of the overlaps
close all

% define the fontsize
fontsize = 15;
% normalize the matrix by the average of the rep distances
% allocate memory for the full control matrix
full_control = zeros(stim_num,stim_num,fish);
% for all the combinations
for combs = 1:number_combinations
    % for all the fish
    for fish = 1:num_fish
        % average per position
        full_control(stim_combinations(combs,1),stim_combinations(combs,2),fish) = ...
            mean(single_matrices{fish,2}([stim_combinations(combs,1),stim_combinations(combs,2)]));
    end
end

% concatenate the single rep results
full_reps = cat(3,single_matrices{:,1});

% allocate memory for the significance
significance_matrix = zeros(stim_num);
% get the significance value for each position, accounting for multiple
% comparisons
for el = 1:stim_num^2
    [x,y] = ind2sub([stim_num,stim_num],el);
    significance_matrix(x,y) = signrank(squeeze(full_reps(x,y,:)-mean(full_control(x,y,:),3))).*number_combinations;
end

figure
imagesc(significance_matrix)
axis square
% set the color scale to the max number of trials per category
if contains(data.name,'p8_S')
    fig_name = 'AF10';
else
    fig_name = data.figure_name;
end
title(fig_name)
set(gca,'CLim',[0, 1])
set(gca,'TickLength',[0 0])
set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',fontsize,...
    'XTickLabelRotation',45)
set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',fontsize)
set(gca,'FontSize',20,'LineWidth',2)
if ~contains(data.name,{'Syn','syn'})
    cba = colorbar;
    set(cba,'TickLength',0,'LineWidth',2)
    ylabel(cba,'Similarity Index')
end
colormap(magma)

% Plot the overlap matrix
close all
% define the fontsize
fontsize = 7;
% average the matrices across animals
norm_overlap = mean(cat(3,single_matrices{:,1}),3);
control_overlap = mean(cat(2,single_matrices{:,2}),2);
% normalize the matrix by the average of the rep distances
% for all the combinations
for combs = 1:number_combinations
    % normalize
    norm_overlap(stim_combinations(combs,1),stim_combinations(combs,2)) = ...
        (norm_overlap(stim_combinations(combs,1),stim_combinations(combs,2))./mean(...
        control_overlap([stim_combinations(combs,1),stim_combinations(combs,2)])));
    % make it symmetric
    norm_overlap(stim_combinations(combs,2),stim_combinations(combs,1)) = ...
        norm_overlap(stim_combinations(combs,1),stim_combinations(combs,2));
end


h = bettercorr(norm_overlap,magma,[],[0 1],significance_matrix);
axis square

% set the color scale to the max number of trials per category
if contains(data.name,'p8_S')
    fig_name = 'AF10';
else
    fig_name = data.figure_name;
end
title(fig_name)
set(gca,'CLim',[0, 1])

set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',fontsize)

if contains(data.name,'p17b')
    rotation = 0;
else
    rotation = 45;
end

set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',fontsize,...
    'XTickLabelRotation',rotation)
colormap(magma)