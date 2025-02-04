%% 0) Compare main and control datasets

% Load the seeds from the downsampled data and the extracted data
% load the data
clearvars
close all
Paths
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).clusters_path;

% define the stimulus time
stim_time = 21:60;

data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
    cone_color_scheme = [0.5 0 0;0 0.5 0;0 0 0.5;0.5 0 0.5];
    stim_labels = {'Red','Green','Blue','UV'};
else
%     color_scheme = distinguishable_colors(6);
    color_scheme = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    cone_color_scheme = [];
    stim_labels = [];
end
% get the number of data sets
num_data = size(data,2);
% Allocate to clusters with GM model and compare

downsampled_traces = data(2).conc_trace;

% get only the stimulus period
downsampled_traces = reshape(downsampled_traces,[],data(2).time_num,data(2).stim_num);
downsampled_traces = reshape(downsampled_traces(:,stim_time,:),[],length(stim_time)*data(2).stim_num);

roi_number = size(downsampled_traces,1);
%define the sPCA parameters to use
bounds_top = 1:length(stim_time):size(downsampled_traces,2);
bounds_bottom = [bounds_top(2:end)-1,size(downsampled_traces,2)];
bounds = [bounds_top;bounds_bottom];
K = ones(1,data(1).stim_num).*4;
t_bins = ones(data(1).stim_num,1).*10;
pca_vec = ones(data(1).stim_num,1).*1;

%define the vector of cluster numbers to try
clu_vec = [];

replicates = 20;

[~,~,~,~,~,f_data] = sPCA_GMM_cluster_Color(downsampled_traces,bounds...
    ,K,t_bins,pca_vec,[],clu_vec,replicates,[],2);

% use the GMM to allocate the downsampled traces
convolved_clusters = cluster(data(1).GMModel,f_data);
% Infer the valid clusters by comparing original and cleaned up idx

% get only the stim time
conc_trace = data(1).conc_trace;
conc_trace = reshape(conc_trace,[],data(1).time_num,data(1).stim_num);
conc_trace = reshape(conc_trace(:,stim_time,:),[],length(stim_time)*data(1).stim_num);

[~,~,~,~,~,f_data] = sPCA_GMM_cluster_Color(conc_trace,bounds...
    ,K,t_bins,pca_vec,[],clu_vec,replicates,[],2);

% get the original idx
original_idx = cluster(data(1).GMModel,f_data);
% get the final idx
final_idx = data(1).idx_clu;
% get a vector with the cluster number
original_clunum = data(1).GMModel.NumComponents;

% allocate memory for the LUT
idx_LUT = zeros(original_clunum,2);
% for all the original clusters
for clu = 1:original_clunum
    idx_LUT(clu,1) = clu;
    idx_LUT(clu,2) = mode(final_idx(original_idx==clu));
end
% Correct the idx for the kernels

for roi = 1:roi_number
    convolved_clusters(roi) = idx_LUT(convolved_clusters(roi)==idx_LUT(:,1),2);
end
%% 1) Compare the clusters allocation between downsampled and the actual data

close all

% assemble a comparison matrix

% allocate memory for the matrix
cluster_matrix = zeros(data(1).clu_num,data(2).clu_num);

% for all the traces
for roi = 1:roi_number

    % get the row coordinate (local clusters)
    x = convolved_clusters(roi);
    % get the Zhou cluster number
    y = data(2).idx_clu(roi);
    % if either is a zero, skip
    if x == 0 || y == 0
        continue
    end
    % add the position of the cluster as a counter
    cluster_matrix(x,y) = cluster_matrix(x,y) + 1;
end

% plot the matrix
imagesc(sortrows(cluster_matrix,'descend'))
cmap = magma(256);
cmap(1,:) = [1 1 1];
colormap(cmap)

set(gcf,'Color','w')
set(gca,'TickLength',[0 0],'FontSize',15)
%% 2) Calculate types and compare

if contains(data(1).name,'p17b')
    
    close all
    % allocate memory to store the matrices
    type_cell = cell(num_data,5);

    % for all of the datasets
    for datas = 1:num_data
        
        % use the gains
        delta_norm = data(datas).delta_norm;
        % get the 10th percentile
        zero_threshold = prctile(abs(delta_norm),10,1);
        % zero the values below a threshold
        delta_norm(abs(delta_norm)<zero_threshold&abs(delta_norm)>0) = 0;
        % turn negatives into -1 and positives into 1
        delta_norm(delta_norm>0) = 1;
        delta_norm(delta_norm<0) = -1;
        
        % quantify the occurrence of each pattern
        [pattern,ia,ic] = unique(delta_norm,'rows');
        
        % get the number of patterns
        pattern_num = length(ia);
        
        % allocate vector for the number
        pattern_counts = zeros(pattern_num,1);
        % count the occurrences
        % for all the patterns
        for pat = 1:pattern_num
            pattern_counts(pat) = sum(ic==pat);
        end
        
        % sort by abundance
        [pattern_counts_sort,sort_idx] = sort(pattern_counts,'descend');
        
        pattern = pattern(sort_idx,:);

        % allocate memory for the colors
        pattern_full = zeros(size(pattern,1),4,3);
        % transform the indexes into colors
        for channel = 1:3
            pattern_full(pattern(:,channel)==1,channel,channel) = 1;
            pattern_full(pattern(:,channel)==0,channel,:) = 1;
            if channel == 1
                pattern_full(pattern(:,4)==1,4,[1 3]) = 1;
                pattern_full(pattern(:,4)==0,4,:) = 1;
            end
 
        end
        
        % store the matrix with the sorted values
        type_cell{datas,1} = pattern;
        type_cell{datas,2} = pattern_counts_sort./sum(pattern_counts_sort);
        type_cell{datas,3} = pattern_full;
        type_cell{datas,4} = sort_idx;
        type_cell{datas,5} = ic;
        
        % eliminate the patterns with only 1 instance
        elim_vector = pattern_counts_sort<2;
        pattern_counts_sort = pattern_counts_sort(~elim_vector);
        pattern_full = pattern_full(~elim_vector,:,:);
        
        figure
        set(gcf,'Color','w')
        subplot(2,1,2)
        image(permute(pattern_full,[2 1 3]))

        set(gca,'YScale','linear','XTick',[],'Visible','off')

        subplot(2,1,1)
        bar((pattern_counts_sort))
        set(gca,'YScale','linear','XTick',[],'Visible','off')
        

        axis tight
             
    end
    %% Plot a single bar plot for type comparison
    close all
    % get the patterns present in the seconda data set that are not in the
    % first
    
    % get the matching indexes for the patterns of 1 wrt 2
    [~,ia,ib] = intersect(type_cell{1,1},type_cell{2,1},'rows','stable');
    
    % define the vectors for plotting
    
    % get the color code
    [~,type_ia,type_ib] = union(type_cell{1,1},type_cell{2,1},'rows','stable');
    color_code = cat(1,type_cell{1,3}(type_ia,:,:),type_cell{2,3}(type_ib,:,:));
    
    % get the number of total types
    total_number = size(color_code,1);
    
    % allocate memory for the counts
    total_counts = zeros(total_number,2);
    
    % put the counts from the original group
    total_counts(1:length(type_cell{1,2}),1) = type_cell{1,2};
    % initialize a counter for the second dataset exclusive types
    counter = 1;
    % for all the types
    for types = 1:total_number
        % check if the type is there
        if sum(ia==types) ~= 0
            % fill in the corresponding count for the second dataset
            total_counts(types,2) = type_cell{2,2}(ib(ia==types));
        elseif types > length(type_cell{1,4}) 
            % find the count of the first exclusive type
            total_counts(types,2) = type_cell{2,2}(type_ib(counter));
            % update the counter
            counter = counter + 1;
            
        end
    end
    
    % plot the results
    figure
    set(gcf,'Color','w')
    subplot(3,1,2)
    image(permute(color_code,[2 1 3]))
    
    set(gca,'XLim',[-0.15 size(total_counts,1)+0.5])
    set(gca,'YScale','linear','XTick',[],'Visible','off')
    
    subplot(3,1,1)
    bar(total_counts(:,1))
    BarPlotBreak(total_counts(:,1),total_counts(2,1)*1.8,total_counts(1,1)*0.9,'Line',0.6,2)
    set(gca,'YScale','linear','XTick',[],'Visible','off')
    axis tight
    % compute the total difference
    disp(sum(abs(total_counts(:,1)-total_counts(:,2))))
    
    % get the coordinates of the max and second to max
    [~,max_idx] = sort(total_counts(:,2),'descend');
    subplot(3,1,3)
    bar(total_counts(:,2))
    BarPlotBreak(total_counts(:,2),total_counts(max_idx(2),2)*1.8,total_counts(max_idx(1),2)*0.9,'Line',0.6,2)
    set(gca,'YScale','linear','XTick',[],'Visible','off','YDir','reverse')
    
    axis tight
    %% Plot the average response (only for delayed dataset)
    
    % select the original p17b datasets plus the delayed one
    if length(data) == 3
        close all
        
        figure
        
        % for all the datasets
        for datas = 1:3
            % find the cluster with the max number of traces
            [~,max_traces] = max(data(datas).clu_number);
            % get the cluster traces
            cluster_traces = normr_1(data(datas).conc_trace(data(datas).idx_clu==max_traces,21:60),0);
            % calculate the average response
            average_response = mean(cluster_traces,1);
            % get the standard deviation
            std_response = std(cluster_traces,0,1);
            % get the x axis
            x = 1:size(average_response,2);
            % plot
            plot(x,average_response)
            hold on
        end
        
    end
    
end