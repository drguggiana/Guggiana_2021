%% Clean up and load 
clearvars
close all

Paths
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;

data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
else
    color_scheme = [];
end

% define the dataset colors (assuming only AF10 vs Tectum)
dataset_colors = paths.afOT_colors;

% define the colormap
cmap = magma;
% Get the region filtering index

% get the number of dataset
num_data = size(data,2);
% allocate memory for the index
index_cell = cell(num_data,1);
fig_name = cell(num_data,1);
% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        region_list = {'AF10'};
        fig_name{datas} = 'AF10';
    else
        region_list = {'R-TcN','R-TcP'};
        fig_name{datas} = 'Tectum';
    end
    % load the anatomy info
    anatomy_info = data(datas).anatomy_info(:,1);
    
    % separate the traces by region
    [region_cell,~] = region_split(data(datas).single_reps,...
        anatomy_info,data(datas).name,region_combination,region_list);
    % rewrite the index vector
    index_cell{datas} = region_cell{3}==1;
    
end
% get the number fo stimuli
stim_num = data(1).stim_num;
% get the number of time bins
time_num = data(1).time_num;
%define the stim labels based on the paradigm

% get the dataset name
stim_name = data(1).name;
%scan for the p17b
if contains(stim_name, 'p17b')
    %if it's p17b
    stim_labels = {'R','G','B','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'R CK','UV CK','R GR','UV GR','R FL','UV FL'};
end
%% 1) Calculate correlation matrices for each data set

close all

% allocate memory to save the labels
correlations = zeros(stim_num,stim_num, num_data);

%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
        %extract only the stim period
        resh_trace = resh_trace(:,21:60,:);

    end
    
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);
    % allocate memory for the output
    corr_perfish = zeros(stim_num,stim_num,num_fish);
    %get the number of time points per stimulus
    t_perstim = size(resh_trace,2);
    % for all the fish
    for fish = 1:num_fish
        % get only the traces from this fish
        resh_fish = resh_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:,:);
        % get the number of traces for this fish
        fish_trace_num = size(resh_fish,1);
        %now reshape again for calculating the correlation across traces for
        %each stimulus
        corr_trace = reshape(resh_fish,fish_trace_num*t_perstim,stim_num);

        %calculate and plot an allvall corr matrix
        corr_perfish(:,:,fish) = corr(corr_trace);
    end
    % calculate the average correlation
    rho = 1-abs(mean(corr_perfish,3));
    
    % plot the correlation matrix
    h = bettercorr(rho,cmap);
    
    
    set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',7,...
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',7)
    set(gca,'CLim',[0,1])
    set(gca,'TickLength',[0 0],'LineWidth',0.05)
    axis square    
    colormap(magma)
    % save the matrix for plotting
    correlations(:,:,datas) = rho;
        
end
%% 2) Calculate correlation in the time dimension

close all

% define the period to take
corr_period = 10:70;
% get the number of points
t_perstim = length(corr_period);
% allocate memory to save the matrices
correlations = zeros(t_perstim,t_perstim,stim_num,num_data);
%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
        %extract only the stim period
        resh_trace = resh_trace(:,corr_period,:);
    end
    
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);

    % allocate memory for the output
    corr_perfish = zeros(t_perstim,t_perstim,num_fish);


    %now reshape again for calculating the correlation across traces for
    %each stimulus
    
    % for all the stimuli
    for stim = 1:stim_num
        figure
        % get the traces for this stim
        stim_trace = resh_trace(:,:,stim);
        % for all the fish
        for fish = 1:num_fish
            % get only the traces from this fish
            corr_trace = stim_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:);
            %calculate and plot an allvall corr matrix
            corr_perfish(:,:,fish) = corr(corr_trace);
        end
        
        %calculate and plot an allvall corr matrix
        rho = mean(corr_perfish,3);
        % store the matrix for later use
        correlations(:,:,stim,datas) = rho;
        imagesc(rho)

        title(stim_labels{stim},'color',color_scheme(stim,:))
        set(gca,'XTick',[],'YTick',[])
        set(gca,'TickLength',[0 0],'LineWidth',0.05)

        axis square
        colormap(cmap)

    end
        
end
%% 3) Plot averages
close all
traces = figure;

% allocate memory to store the cumulative deltas
delta_matrix = zeros(stim_num,1);
% define the offset
offset = 0.5;
for stim =1:stim_num

    % get the data
    matrix1 = squeeze(correlations(:,:,stim,1));
    matrix2 = squeeze(correlations(:,:,stim,2));
    
    % blank the diagonal
    matrix1(eye(size(matrix1))==1) = NaN;
    matrix2(eye(size(matrix1))==1) = NaN;
    
    mean1 = nanmean(matrix1,1);
    mean2 = nanmean(matrix2,1);
    
    std1 = nanstd(matrix1,0,1)./sqrt(size(matrix1,1)-1);
    std2 = nanstd(matrix2,0,1)./sqrt(size(matrix1,1)-1);
    
    % get the x_range
    x_range = ((1:size(matrix1,1))-10)./data(datas).framerate;
    
    % select the target figure
    figure(traces)
    tectum_handles = shadedErrorBar(x_range,mean1+(stim_num-(stim)*offset),std1,{'-','Color',color_scheme(stim,:)});
    hold on
    af10_handles = shadedErrorBar(x_range,mean2+(stim_num-(stim)*offset),std2,{'--','Color',color_scheme(stim,:)});
    % color the line black
    set(tectum_handles.mainLine,'Color',[0 0 0],'LineWidth',1);
    set(af10_handles.mainLine,'Color',[0 0 0],'LineWidth',1);

end

figure(traces)

set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',15)
set(gca,'YTick',[])
xlabel('Time (s)')
ylabel('Relative Corr.')
box off
axis tight
plot([0 0],get(gca,'YLim'),'k--','LineWidth',1)
%% 4) Calculate correlation over time for p8

% if it's not p8, skip
if contains(data(1).name,'p8')
    close all

    % define the marker size
    marker_size = 2;
    markers = {'s','o'};
    % marker shape left
    shape_left_x = {[-0.1 0 0],[-0.1 0 0 -0.1]};
    shape_left_y = {[-0.01 -0.01 0.01],[-0.01 -0.01 0.01 0.01]};
    % marker shape right
    shape_right_x = {[0.1 0 0],[0.1 0 0 0.1]};
    shape_right_y = {[-0.01 -0.01 0.01],[-0.01 -0.01 0.01 0.01]};
    %define the sets of time regions to correlate
    time_corr = (1:50)';

    %get the time axis
    %get the number of time points
    timep_num = size(time_corr,1);
    %allocate memory for the axis
    timep_axis = zeros(timep_num,1);
    %for all the time points
    for timep = 1:timep_num
        timep_axis(timep) = mean(time_corr(timep,:));
    end
    %get the number of time regions
    num_times = size(time_corr,1);

    %allocate memory to store the correlations
    tcorr_mat = zeros(num_data,num_times,stim_num,stim_num);
    tcorr_mat_sem = zeros(num_data,num_times,stim_num,stim_num);
    % generate a single figure for the subplots
    figure

    %for both data sets
    for datas = 1:num_data

        %if a normal data set
        if datas <= num_data
            %load the clusters
            conc_trace = data(datas).conc_trace;
            %get the number of traces
            trace_num = size(conc_trace,1);
            %reshape the matrix to calculate correlations between stimuli
            resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
            %extract only the stim period
            resh_trace = resh_trace(:,11:60,:);
        end
        % get the trace fish of origin
        fish_ori = data(datas).fish_ori;
        % get the number of fish
        num_fish = size(unique(fish_ori(:,1)),1);

        %for all the times
        for times = 1:num_times

            % allocate memory for the output
            corr_perfish = zeros(stim_num,stim_num,num_fish);

            %now reshape again for calculating the correlation across traces for
            %each stimulus
            corr_trace = squeeze(mean(resh_trace(:,time_corr(times,:),:),2));
            % for all the fish
            for fish = 1:num_fish
                % get only the traces from this fish
                resh_fish = corr_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:,:);

                %calculate and plot an allvall corr matrix
                corr_perfish(:,:,fish) = corr(resh_fish);
            end
            % save the average and sem corr matrix
            tcorr_mat(datas,times,:,:) = mean(corr_perfish,3);
            tcorr_mat_sem(datas,times,:,:) = std(corr_perfish,0,3)./sqrt(num_fish);
        end
        
        % define the combinations manually
        comb_vec = [1 2;3 4;5 6];
        colors_reduv = [1 0 0;1 0 1];
        time_axis = (1:size(tcorr_mat,2));
        %get the number of combinations
        comb_num = size(comb_vec,1);

        %allocate memory for the legend
        legend_cell = cell(comb_num,1);
        %for all the combs
        for combs = 1:comb_num
            dat1 = squeeze(tcorr_mat(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
            sem1 = squeeze(tcorr_mat_sem(datas,:,comb_vec(combs,1),comb_vec(combs,2)));

            % select the correct subplot
            subplot(3,1,combs)
            
            hold('on')

            errorbar(time_axis,dat1,sem1,strcat('-o'),'MarkerSize',3,...
                'CapSize',2,'MarkerEdgeColor',dataset_colors(datas,:),...
                'Color',dataset_colors(datas,:))
            %assemble the legend
            legend_cell{combs} = strcat(stim_labels{comb_vec(combs,1)},'_',stim_labels{comb_vec(combs,2)});
            set(gca,'FontSize',7)
            set(gca,'XLim',[0 size(time_axis,2)],'YLim',[-0.2 0.7],'YTick',[0 0.3 0.6])
            
            switch combs
                case 3
                    xlabel('Time (s)','FontSize',7)
                    ylabel('Flash')
                case 2
                    set(gca,'XTick',[])
                    ylabel('Grating')
                case 1
                    set(gca,'XTick',[])
                    ylabel('Checker')
            end

            set(gca,'TickLength',[0 0])
            plot([10 10],get(gca,'YLim'),'k-')
            box off
        end

    end
end
%% 5) Calculate decorrelation matrix and histograms

close all

% get the stim_number
stim_num = data(1).stim_num;

% create the figure and get the handle
[rho_cell,combo_list,combo_number] = matrix_subplot(data,index_cell);

% allocate memory for a matrix
result_matrix = zeros(stim_num);
% and for the tests
test_matrix = zeros(combo_number,1);

% for all the combinations
for combo = 1:combo_number

    % load the data
    data1 = 1-abs(rho_cell{combo_list(combo,1),combo_list(combo,2),1});
    data2 = 1-abs(rho_cell{combo_list(combo,1),combo_list(combo,2),2});
    
    result_matrix(combo_list(combo,1),combo_list(combo,2)) = median(data1(:));
    result_matrix(combo_list(combo,2),combo_list(combo,1)) = median(data2(:));
    
    % test, adjusting for multiple comparisons
    test_matrix(combo) = ranksum(abs(data1(:)),abs(data2(:))).*combo_number;

end

% plot the matrix
figure
imagesc(result_matrix)
axis square
set(gca,'TickLength',[0 0],'LineWidth',0.05)
set(gca,'XTick',1:stim_num,'YTick',1:stim_num,'XTickLabels',...
    stim_labels,'YTickLabels',stim_labels,'XTickLabelRotation',45) 
set(gca,'FontSize',15)

cmap = magma;
cmap(1,:) = [1 1 1];
colormap(cmap)
%% 6) Correlate the responses across fish

close all
% define the font size
fontsize = 20;
%for both data sets
for datas = 1:num_data
    figure
    %load the clusters
    conc_trace = data(datas).conc_trace;
    
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);
    % get the number of clusters
    clu_num = data(datas).clu_num;
    % allocate memory for the averages
    fish_average = zeros(clu_num,size(conc_trace,2),num_fish);
    fish_cluster_traces = zeros(clu_num,num_fish);
    % calculate the per-fish cluster averages
    for fish = 1:num_fish
        % get the idx_clu for this fish
        idx_fish = data(datas).idx_clu(fish_ori(:,1)==fish);
        % get the traces for this fish
        fish_trace = conc_trace(fish_ori(:,1)==fish,:);
        % for all the clusters
        for clu = 1:clu_num
            fish_average(clu,:,fish) = nanmean(fish_trace(idx_fish==clu,:),1);
            fish_cluster_traces(clu,fish) = sum(idx_fish==clu);
        end
    end
    % get all the combinations of fish
    fish_comb = nchoosek(1:num_fish,2);
    number_combos = length(fish_comb);
    % allocate memory for a correlation matrix
    fish_correlation = zeros(num_fish);
    % for all the fish combos
    for combo = 1:number_combos

        % get the cluster averages
        fish1 = fish_average(:,:,fish_comb(combo,1));
        fish2 = fish_average(:,:,fish_comb(combo,2));
        % get the corelation matrix
        corr_matrix = corr(fish1',fish2');
        % take the off diagonal entries (since clusters are matched)
        corr_matrix = diag(corr_matrix);
        % average
        fish_correlation(fish_comb(combo,1),fish_comb(combo,2)) = ...
            nanmean(corr_matrix(:));
        
    end
    % plot the matrix
    imagesc(fish_correlation)
    set(gca,'TickLength',[0 0],'LineWidth',2,'XTick',1:num_fish,'YTick',1:num_fish)
    set(gca,'CLim',[0 1])
    xlabel('Fish')
    ylabel('Fish')
    
    axis square
    axis tight
    colormap(cmap)

    % Plot the number of traces of every cluster in every fish
    figure
    imagesc(log(fish_cluster_traces./max(fish_cluster_traces,[],1)))

    set(gca,'TickLength',[0 0],'LineWidth',2,'XTick',1:num_fish,'YTick',1:3:size(fish_cluster_traces,1))
    xlabel('Fish')
    ylabel('Cluster')
    
    axis square
    axis tight
    colormap(cmap)
    
end