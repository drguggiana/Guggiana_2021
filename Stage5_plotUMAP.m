%% 0) Embed the data with UMAP
clearvars
close all

Paths
addpath(genpath(paths(1).main_path))
main_path = paths(1).clusters_path;

% reset the rng
rng(1);

data = load_clusters(main_path);

% get the number of dataset
num_data = size(data,2);
% allocate memory for the index
index_cell = cell(num_data,1);
% allocate memory for the region list
region_list = cell(num_data,1);
% define the region set to use (1 tc vs rgc, 2 all vs all)
region_set = 1;
hist_colors = zeros(num_data,3);
% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        if region_set == 1
            region_list{datas} = {'AF10'};
            hist_colors(datas,:) = paths.afOT_colors(datas,:);
        elseif region_set == 2
            region_list{datas} = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
        end
    else
        if region_set == 1
            region_list{datas} = {'R-TcN','R-TcP'};
            hist_colors(datas,:) = paths.afOT_colors(datas,:);
        elseif region_set == 2
            region_list{datas} = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
        end
    end
    % load the anatomy info
    anatomy_info = data(datas).anatomy_info(:,1);
    
    % separate the traces by region
    [region_cell,~] = region_split(data(datas).single_reps,...
        anatomy_info,data(datas).name,region_combination,region_list{datas});
    % rewrite the index vector
    index_cell{datas} = region_cell{3}==1;
    
end

% define the labels for the gain histogram
gain_labels = {'NR','Red','Green','Blue','UV','2-chrom','3-chrom','4-chrom'};

% Allocate memory to store the embeddings
UMAP_cell = cell(num_data,1);
% for all the datasets selected
for datas = 1:length(data)
    close all
    %scan for the p17b
    if contains(data(datas).name,'p17b')
        %if it's p17b
        stim_labels = {'Red','Green','Blue','UV'};
        %define the plot colors
        plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];
        % set gain to on
        gains = 1;
        % get the gains
        delta_norm = data(datas).delta_norm(index_cell{datas}==1,:);
    else %if it's p6p8 instead
        %define the stim labels (for p6p8 data)
        stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
        %define the plot colors
%         plot_col = [1 0 0;0 0 1;0.8 0 0;0 0 0.8;0.4 0 0;0 0 0.4];
        plot_col = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
        gains = 0;
    end
    %% Extract the PCA features from the data
    
    % get the PC info from the structure
    pcs = data(datas).pcs(:,1);
    % get the data
    conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);

    % allocate memory for the vector
    pc_matrix = zeros(size(conc_trace,1),size(pcs{1},2),data(datas).stim_num);
    % define the stimulus time
    stim_time = 21:60;
    % get the data(datas) and reshape to separate stim and time
    stim_data = reshape(conc_trace,[],data(datas).time_num,data(datas).stim_num);
    % for all the stimuli
    for stim = 1:data(datas).stim_num
        pc_matrix(:,:,stim) = stim_data(:,stim_time,stim)*pcs{stim};
    end
    % concatenate and column normalize the matrix for use
    pc_matrix = normr_1(reshape(pc_matrix,[],size(pcs{1},2)*data(datas).stim_num),2);
    %% Use UMAP to embed the gains
    
    [reduced_data, umap] = run_umap(pc_matrix, 'n_neighbors', 10, 'min_dist', 0.1);
%     [reduced_data, umap] = run_umap(cat(2,pc_matrix,normr_1(data(datas).delta_norm,2)), 'n_neighbors', 10, 'min_dist', 0.1);
    % store the embedding
    UMAP_cell{datas} = reduced_data;
end
%% 1) Plot the results based on stimulus
close all
histo = figure;

% define whether to use gains or average levels
metric = 'average';
% allocate memory to store the max values for stats
max_cell = cell(num_data,1);
for datas = num_data:-1:1
    
    % load the embedding
    reduced_data = UMAP_cell{datas};
    num_points = size(reduced_data,1);
    % get the number of traces, time and stimuli
    conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);
    
    trace_num = size(conc_trace,1);
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    
    switch metric
        case 'average'
            % get the reshaped activity
            average_levels = reshape(conc_trace,trace_num,time_num,stim_num);
            % take only the stimulation time
            average_levels = average_levels(:,stim_time,:);
            % take the absolute average
            average_levels = squeeze(mean(abs(average_levels),2));
        case 'gains'
            % use the actual gains
            average_levels = abs(data(datas).delta_norm(index_cell{datas}==1,:));
    end

    if gains

        % use the top percentile of the gain
        perc_values = prctile(average_levels,75,1);
        % allocate memory for the max_values
        max_values = zeros(size(average_levels,1),1);
        temp_values = zeros(size(average_levels));
        % for all colors
        for stim = 1:size(average_levels,2)
            % get the ROIs passing threshold in each color
            temp_values(average_levels(:,stim)>perc_values(stim),stim) = 1;
            
            max_values(average_levels(:,stim)>perc_values(stim)) = stim;
        end
        
        % get the ROIs with more than 1 selectivity
        opponents = sum(temp_values,2)+3;
        opponents(opponents<5) = 0;
        % join it with the max_values
        max_values(opponents~=0) = opponents(opponents~=0);
        max_values = max_values + 1;
        % store the value for stats
        max_cell{datas} = max_values;


        if datas == 1
            figure(histo)
            
            % allocate memory for the bin counts
            bin_counts = cell(2,1);
            % get the histcounts
            for i = 1:2
                bin_counts{i} = histcounts(max_cell{i},'Normalization','probability');
            end
            bars = bar(vertcat(bin_counts{:})');
            tint_colors = tint_colormap(hist_colors,0);
            set(bars(1),'FaceColor',tint_colors(2,:))
            set(bars(2),'FaceColor',tint_colors(1,:))
            
            box off
            set(gca,'XTick',1:8,'XTickLabels',gain_labels,'XTickLabelRotation',45)
            set(gca,'FontSize',7,'LineWidth',1,'TickLength',[0 0])

        end
        
        figure
        
        cmap = [0.8 0.8 0.8;1 0 0;0 1 0;0 0 1;1 0 1;0 0 0;0 0 0;0 0 0];
        gscatter(reduced_data(:,1),reduced_data(:,2),max_values,cmap,'.',10)

        legend('off')
        set(gca,'TickLength',[0 0],'visible','off','LineWidth',2)
        set(gcf,'Color','w')
        title(strjoin({'UMAP',data(datas).name,'Combined'},'_'),'Interpreter','None')

    else

        
        figure
        [~,max_values] = max(average_levels,[],2);

        cmap = [1 0 0;1 0 1;1 0 0;1 0 1;0.5 0 0;0.5 0 0.5];
        emap = [1 0 0;1 0 1;0 0 0;0 0 0;0.5 0 0;0.5 0 0.5];
        markers = {'o','o','s','s','d','d'};
        for stim = 1:stim_num
            
            scatter(reduced_data(max_values==stim,1),reduced_data(max_values==stim,2),...
                [],markers{stim},'MarkerEdgeColor',emap(stim,:),'MarkerFaceColor',cmap(stim,:))
            hold on
        end
        legend({'Red CK', 'UV CK', 'Red GR', 'UV GR', 'Red FL', 'UV FL'},'Location','northwest')
        colormap(cmap)
        colorbar('Ticklabels',stim_labels,'TickLabelInterpreter','None','Ticks',linspace(0.08,0.92,6))
        title(strjoin({data(datas).name,'Combined'},'_'),'Interpreter','None')
        set(gca,'TickLength',[0 0],'visible','off')

    end
end
% Run stats on the histogram of types

% allocate memory for the counts for both data sets
full_counts = cell(num_data,1);
% for both data sets
for datas = 1:num_data
    % get the fish info
    fish_info = data(datas).fish_ori(index_cell{datas}==1,1);
    % get the number of fish
    num_fish = length(unique(fish_info));
    
    % get the types of bins and their number
    bin_types = unique(max_cell{datas});
    type_number = length(bin_types);
    
    % allocate memory for the count
    count_per_fish = zeros(num_fish,type_number);
    
    % for the types of bins
    for bintype = 1:type_number
        % for all the fish
        for fish = 1:num_fish
            % get the current max values
            current_values = max_cell{datas}(fish_info==fish);
            count_per_fish(fish,bintype) = sum(current_values==bin_types(bintype));
        end
    end
    % store the results
    full_counts{datas} = count_per_fish;
end

% allocate memory for the test results
test_result = zeros(type_number,1);
% run wilcoxon's sign rank for each type
for bintype = 1:type_number
    test_result(bintype) = signrank(full_counts{1}(:,bintype),full_counts{2}(:,bintype));
end
%% 2) Plot the p8 UMAPs separately by stimulus modality
close all
if contains(data(datas).name,'p8')
    % define the labels of the patterns
    pattern_label = {'Checker','Gratings','Flash'};
    % define the dataset names
    fig_name = {'Tectum','AF10'};
    % allocate memory to store the color raw vectors
    color_cell = cell(num_data,stim_num/2);
    % for all datasets
    for datas = 1:num_data

        % load the embedding
        reduced_data = UMAP_cell{datas};
        num_points = size(reduced_data,1);
        % close all
        % get the number of traces, time and stimuli
        conc_trace = data(datas).conc_trace;
        trace_num = size(conc_trace,1);
        time_num = data(datas).time_num;
        stim_num = data(datas).stim_num;

        % get the reshaped activity
        average_levels2 = reshape(conc_trace,trace_num,time_num,stim_num);
        % take only the stimulation time
        average_levels2 = average_levels2(:,stim_time,:);
        % take the absolute average
        average_levels2 = squeeze(max(abs(average_levels2),[],2));

        % define the colormap
        cmap = [0.8 0.8 0.8;1 0 0;1 0 1;0 0 0];
        % for all stimulus types
        for stim = 1:2:5
            figure
            
            % copy the corresponding average_levels
            average_levels = average_levels2(:,stim:stim+1);
            
            perc_values = prctile(average_levels,75,1);
            % allocate memory for the max_values
            max_values = zeros(size(average_levels,1),1);
            temp_values = zeros(size(average_levels));
            % for all colors
            for stim2 = 1:2
                % get the ROIs passing threshold in each color
                temp_values(average_levels(:,stim2)>perc_values(stim2),stim2) = 1;
                
                max_values(average_levels(:,stim2)>perc_values(stim2)) = stim2;
            end
            
            % get the ROIs with more than 1 selectivity
            opponents = sum(temp_values,2)+1;
            opponents(opponents<3) = 0;
            % join it with the max_values
            max_values(opponents~=0) = opponents(opponents~=0);
            color_raw = max_values + 1;

            % save the color_raw vector for quantification
            color_cell{datas,(stim+1)/2} = color_raw;
            % plot the scatter
            s = gscatter(reduced_data(:,1),reduced_data(:,2),color_raw,cmap,'.',...
                1,'off');

            colormap(cmap)
            set(gcf,'Color','w')
            title(strjoin({data(datas).name,'Combined'},'_'),'Interpreter','None')
            set(gca,'TickLength',[0 0],'visible','off')
        end
        %% Produce a venn diagram
        close all
        figure
        
        fontsize = 15;
        % turn the color matrix into binary, excluding mixed selectivity
        % cells
        binary_matrix = color_all(~any(color_all>3,2),:);
        binary_matrix = binary_matrix>1;
        % get the number of multicolor cells
        multi_number = size(color_all,1)-size(binary_matrix,1);
        
        % get the unique patterns
        [unique_seq,ia,ic] = unique(binary_matrix,'rows');
        % allocate memory for the amounts of ROIs
        roi_counts = zeros(size(unique_seq,1),1);
        % get the counts per pattern
        for pattern = 1:size(unique_seq)
            roi_counts(pattern) = sum(ic==pattern);
        end
        % reorder the pattern
        roi_counts = roi_counts([1 2 3 5 4 6 7 8]);
        unique_seq = unique_seq([1 2 3 5 4 6 7 8],:);
        % draw the venn diagram
        
        colors = magma(3);
        [H,S] = venn(roi_counts(2:end),'FaceColor',{colors(1,:), colors(2,:), colors(3,:)});
        hold on
        
        
        % get the normalized numbers as percentages
        roi_percentages = round(roi_counts*100./sum(roi_counts));
        
        % place text on the circles
        % for all the areas
        for areas = 1:7
            text(S.ZoneCentroid(areas,1),S.ZoneCentroid(areas,2),...
                strcat(num2str(roi_percentages(areas+1)),'%'),...
                'FontSize',7,'HorizontalAlignment','center')
        end
        
        % add the NR and Multi text
        text(min(S.ZoneCentroid(:,1))-3.5,max(S.ZoneCentroid(:,2)),strcat('NR:',num2str(roi_percentages(1)),'%'),...
             'FontSize',7,'HorizontalAlignment','center')
        
        set(gca,'Visible','off')
        axis equal
    end
    %% Spectral types
    close all
    % allocate memory for storing counts for stats
    type_cell = cell(num_data,5);
    % define the color map for the types
    type_colors = [1 1 1;1 0 0;1 0 1;0.5 0 0.5];
    for datas = 1:num_data
        
        color_all = horzcat(color_cell{datas,:});
        
        % remove the non-responsive
        color_responsive = color_all(~(sum(color_all,2)==3),:);
        % quantify the occurrence of each pattern
        [pattern,ia,ic] = unique(color_responsive,'rows');
        
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
        pattern_full = zeros(size(pattern,1),3,3);
        % transform the indexes into colors
        for modality = 1:3
            for stim = 1:4
                pattern_full(pattern(:,modality)==stim,modality,1) = type_colors(stim,1);
                pattern_full(pattern(:,modality)==stim,modality,2) = type_colors(stim,2);
                pattern_full(pattern(:,modality)==stim,modality,3) = type_colors(stim,3);
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
        bar((pattern_counts_sort),'FaceColor',[0.8 0.8 0.8])
        set(gca,'YScale','linear','XTick',[],'Visible','off')
        

        axis tight
       
    end
    %% Compare the histograms statistically
    
    % get the histograms
    h1 = histcounts(type_cell{1,5},64,'Normalization','Probability');
    h2 = histcounts(type_cell{2,5},64,'Normalization','Probability');
    pdist2_piotr(h1,h2)
    
end