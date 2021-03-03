%% Cluster groups of fish together
%% clean up
clearvars
close all force
matlabrc
Paths
addpath(genpath(paths(1).main_path))
%% Load the files and define paths

% reset the rng to keep results reproducible
rng(3);

%get the folder where the image files are
tar_path_all = uipickfiles('FilterSpec',...
    paths(1).rois_path);

%get the number of experiments selected
num_exp = length(tar_path_all);

%define the list of labels to sort the files
label_list = {'_thres.mat'};

%get the program used
[~,name_whole,~] = fileparts(tar_path_all{1});
name_parts = strsplit(name_whole,'_');

%get the number of each type of data (the round is to avoid the nonscalar
%warning for the for loop)
num_data = round(num_exp./length(label_list));

%allocate memory for the different types of files
name_cell = cell(num_data,length(label_list));

%for the types of files
for f_type = 1:length(label_list)
    %get the coordinates of the file names
    name_map = ~cellfun(@isempty,strfind(tar_path_all,label_list{f_type}));
    %store them in the corresponding layer of the name cell
    name_cell(:,f_type) = tar_path_all(name_map);
end
%% Load the traces (prep for clustering)

% for all the files
for files = 1:num_data
    
    %% Define/load constants
    % define the name to save the files
    [~,ori_name] = fileparts(name_cell{files,1});
    ori_name = ori_name(1:end-6);

    %get the number of time points
    time_num = load(name_cell{files,1},'time_num');
    time_num = time_num.time_num;
    %get the number of stimuli
    stim_num2 = load(name_cell{files,1},'stim_num2');
    stim_num2 = stim_num2.stim_num2;
    
    % define the stimulus period (in frames) within a trial
    stim_time = 21:60;
    % load the color info
    col_out = load(name_cell{files,1},'col_out');
    col_out = col_out.col_out;

    %load and concatenate all the traces

    %load the conc_trace variable into the trace_all cell
    conc_trace = load(name_cell{files},'conc_trace');
    %extract the matrix of values from the struct
    conc_trace = conc_trace.conc_trace;

    %load the fish origin information
    fish_ori = load(name_cell{files},'fish_ori');
    fish_ori = fish_ori.fish_ori;

    % load the anatomy info
    anatomy_info = load(name_cell{files},'cat_anatomy_all');
    anatomy_info = anatomy_info.cat_anatomy_all;
    
    % load the single reps
    single_reps = load(name_cell{files},'cat_reps');
    single_reps = single_reps.cat_reps;

    %get the number of fish
    fish_num = length(unique(fish_ori(:,1)));

    %get the number of traces
    trace_num = size(conc_trace,1);
    %% Exclude the outsides of the trial (pre and post stim)

    %reshape the matrix for conversion
    conc_trace2 = reshape(conc_trace,trace_num,time_num,[]);
    %eliminate the undesired time fragments
    conc_trace2 = conc_trace2(:,stim_time,:);
    %rewrite the time_num and conc_trace variables
    time_num2 = size(conc_trace2,2);
    conc_trace2 = reshape(conc_trace2,trace_num,[]);
    %% assemble a structure with all of the data
    
    % create the structure
    main_str = struct([]);
    
    % fill in basic values
    main_str(1).name = ori_name;
    main_str(1).conc_trace = conc_trace;
    main_str(1).single_reps = single_reps;
    main_str(1).col_out = col_out;
    main_str(1).stim_num = stim_num2;
    main_str(1).time_num = time_num;
    %% Extract the gains
    if contains(ori_name,'p17b')
        % get the gains
        [delta_norm, qual_res, cross_res] = gain_analysis(main_str,stim_time,paths.param_path);
        main_str(1).delta_norm = delta_norm;
        main_str(1).qual_res = qual_res;
        main_str(1).cross_res = cross_res;
    else
        main_str(1).delta_norm = [];
        main_str(1).qual_res = [];
        main_str(1).cross_res = [];
    end    
    %% Cluster the traces
    %define the sPCA parameters to use
    bounds_top = 1:time_num2:size(conc_trace2,2);
    bounds_bottom = [bounds_top(2:end)-1,size(conc_trace2,2)];
    bounds = [bounds_top;bounds_bottom];
    K = ones(1,stim_num2).*4;
    t_bins = ones(stim_num2,1).*10;
    pca_vec = ones(stim_num2,1).*1;

    %define the vector of cluster numbers to try
    clu_vec = [5 10 20 30 50 70];
    
    replicates = 20;
    
    [idx_clu,GMModel,clu_num,pcs,bic_vec] = sPCA_GMM_cluster_Color(conc_trace2,bounds...
        ,K,t_bins,pca_vec,[],clu_vec,replicates);

    %save the original cluster number for the surrogate analysis
    ori_clu = clu_num;
    %% Use correlation to merge clusters that are too similar
    % define the correlation threshold
    corr_threshold = 0.8;
    close all
    %allocate memory for the averages
    clu_ave = zeros(clu_num,size(conc_trace,2));
    %and for the trace number
    clu_number = zeros(clu_num,1);
    %for all the clusters
    for clu = 1:clu_num
        %calculate the cluster average
        clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
        %and store the number of traces going into each average
        clu_number(clu) = sum(idx_clu==clu);
    end
    
    % calculate the correlation matrix between clusters
    [cluster_corr,pval] = corr(clu_ave');
    % keep only the significant correlations
    cluster_corr(pval>0.05) = NaN;
    
    figure
    imagesc(cluster_corr)
    figure
    temp_corr = tril(cluster_corr,1);
    histogram(temp_corr(:),'Normalization','cdf')

    % get the matrix with correlation values
    corr_submatrix = triu(cluster_corr,1);
    % allocate memory for the new idx
    new_idx = idx_clu;
    % go through the clusters
    for clu = 1:clu_num
        % check if this cluster has already been eliminated
        if all(isnan(corr_submatrix(clu,:)))
            continue
        end
        % get the corresponding row
        curr_corr = corr_submatrix(clu,:);
        % nan the matching cluster
        curr_corr(clu) = NaN;
        % get the clusters that match with this one above threshold
        matching_clusters = find(curr_corr>corr_threshold);
        % for all the matching clusters
        for match = 1:length(matching_clusters)
            new_idx(new_idx==matching_clusters(match)) = clu;
        end
    end
    % renumber the clusters consecutively
    % get the new cluster numbers
    new_numbers = unique(new_idx);
    % get the new clu_num
    new_clu_num = length(new_numbers);
    % renumber
    for clu = 1:new_clu_num
        new_idx(new_idx==new_numbers(clu)) = clu;
    end
    % redefine the cluster number and idx
    idx_clu = new_idx;
    clu_num = new_clu_num;
    %% Exclude clusters with snr below a threshold and with too few traces
    
    % define the thresholds
    num_thres = 10;
    stim_thres = 0;

    %load the snr for the dataset
    snr_mat = load(name_cell{files},'snr_mat');
    snr_mat = snr_mat.snr_mat;
    
    [idx_clu,clu_num] = cluster_snr(snr_mat,clu_num,idx_clu,num_thres,stim_thres);
    %% Renumber clusters based on number of traces
    
    % get the numbers of traces for each cluster
    [idx_count,edges] = histcounts(idx_clu,(1:clu_num+1)-0.5);
    % sort and get the LUT
    [~,idx_sort] = sort(idx_count,'descend');
    % allocate memory for the replacement index
    idx_new = zeros(size(idx_clu));
    % for all the clusters
    for clu = 1:clu_num
        idx_new(idx_clu==idx_sort(clu)) = clu;
    end
    % replace the original idx_clu
    idx_clu = idx_new;
    %% Calculate the cluster average
    
    %allocate memory for the averages
    clu_ave = zeros(clu_num,size(conc_trace,2));
    %and for the trace number
    clu_number = zeros(clu_num,1);
    %for all the clusters
    for clu = 1:clu_num
        %calculate the cluster average
        clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
        %and store the number of traces going into each average
        clu_number(clu) = sum(idx_clu==clu);
    end
    %% Cluster by area and exclude noisy clusters
    
    % if there's no anatomy info, skip
    if ~isempty(anatomy_info)
        
  
        % get the number of regions
        region_set = unique(anatomy_info(:,1));
        number_regions = length(region_set(~isnan(region_set)));

        % allocate memory to store the region clusters
        region_clusters = struct([]);
        % for all the regions
        for regions = 1:number_regions
            close all
            % get the traces for the region
            region_traces = conc_trace2(anatomy_info(:,1)==region_set(regions),:);
            if contains(ori_name,'p17b')
                % also the gains
                region_gains = main_str(1).delta_norm(anatomy_info(:,1)==region_set(regions),:);
            else
                region_gains = [];
            end
                % is there's less than num_thres traces, skip
                if size(region_traces,1) < num_thres || size(region_traces,1) < sum(K)
                    region_clusters(regions).idx_clu = [];
                    region_clusters(regions).clu_num = [];
                    region_clusters(regions).clu_ave = [];
                    region_clusters(regions).trace_num = [];
                    continue
                end
                % cluster them as above
                [idx_clu_region,~,clu_num_region] = sPCA_GMM_cluster_Color(region_traces,bounds...
                ,K,t_bins,pca_vec,[],clu_vec,replicates,region_gains);

                % filter the traces based on the above
                [idx_clu_region,clu_num_region] = cluster_snr(snr_mat(anatomy_info(:,1)==region_set(regions),:),...
                    clu_num_region,idx_clu_region,num_thres,stim_thres);
                
                % if no clusters were left, also skip
                if sum(idx_clu_region) == 0
                    region_clusters(regions).idx_clu = [];
                    region_clusters(regions).clu_num = [];
                    region_clusters(regions).clu_ave = [];
                    region_clusters(regions).trace_num = [];
                    continue
                end

                % calculate the cluster averages
                %allocate memory for the averages
                clu_ave_region = zeros(clu_num_region,size(region_traces,2));
                %and for the trace number
                clu_number_region = zeros(clu_num_region,1);
                %for all the clusters
                for clu = 1:clu_num_region
                    %calculate the cluster average
                    clu_ave_region(clu,:) = mean(region_traces(idx_clu_region==clu,:),1);
                    %and store the number of traces going into each average
                    clu_number_region(clu) = sum(idx_clu_region==clu);
                end

                % save the cluster info
                region_clusters(regions).idx_clu = idx_clu_region;
                region_clusters(regions).clu_num = clu_num_region;
                region_clusters(regions).clu_ave = clu_ave_region;
                region_clusters(regions).trace_num = clu_number_region;
        end 
    else
        region_clusters(regions).idx_clu = [];
        region_clusters(regions).clu_num = [];
        region_clusters(regions).clu_ave = [];
        region_clusters(regions).trace_num = [];
    end
    %% Load the info about seeds
    
    z_seed = load(name_cell{files,1},'cat_z_all');
    z_seed = z_seed.cat_z_all;
    
    xy_seed = load(name_cell{files,1},'cat_seed_all');
    xy_seed = xy_seed.cat_seed_all;
    
    
    ave_stack = load(name_cell{files,1},'cat_stack_all');
    ave_stack = ave_stack.cat_stack_all;
    %% Add the framerate and official name based on the dataset
    switch ori_name
        case 'p8_gc6s'
            framerate = 1./0.736;
            figure_name = 'Tectum';
        case 'p8_SynG6s'
            framerate = 1./0.507;
            figure_name = 'AF10';
        case 'p17b_gc6s'
            framerate = 1./0.952;
            figure_name = 'RAs';
        case 'p17b_syngc6s'
            framerate = 1./0.952;
            figure_name = 'RGCs';
        case 'p17b_h2b6s'
            framerate = 1./0.952;
            figure_name = 'Tectum, Nuclear labelled';
        case 'p17bdownsample_gc6s'
            framerate = 1/0.952;
            figure_name = 'Downsampled RAs';
        otherwise
            if contains(ori_name,'conv')&&contains(ori_name,'p17b')
                framerate = 1/0.952;
                figure_name = 'Delayed RGCs';
            else
                framerate = 0;
                figure_name = 'Other';
            end
                
            
    end
    %% Assemble the structure with the data

    % continue filling the structure
    main_str(1).fish_ori = fish_ori;
    main_str(1).anatomy_info = anatomy_info;
    main_str(1).pcs = pcs;
    main_str(1).GMModel = GMModel;
    main_str(1).idx_clu = idx_clu;
    main_str(1).bounds = bounds;
    main_str(1).K = K;
    main_str(1).t_bins = t_bins;
    main_str(1).pca_vec = pca_vec;
    main_str(1).bic_vec = bic_vec;
    main_str(1).clu_num = clu_num;
    main_str(1).clu_number = clu_number;
    main_str(1).clu_ave = clu_ave;
    main_str(1).z_seed = z_seed;
    main_str(1).xy_seed = xy_seed;
    main_str(1).ave_stack = ave_stack;
    main_str(1).snr_mat = snr_mat;
    main_str(1).region_clusters = region_clusters;
    main_str(1).framerate = framerate;
    main_str(1).figure_name = figure_name;
    %% Save the structure
    
    %define the save path
    save_path = paths(1).clusters_path;

    % assemble the file name
    save_clu = strcat(ori_name,'_clusters.mat');
    
    % save the structure
    save(fullfile(save_path,save_clu),'main_str')
end