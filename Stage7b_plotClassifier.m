%% 0) Plot classifier results

clearvars
close all
Paths
addpath(genpath(paths(1).main_path))
classifier_path = paths(1).classifier_path;
data_struct = load_clusters(classifier_path);

% get the file names
file_names = {data_struct.name}';
classpcolor = cat(1,data_struct.classpcolor);
shuff_label = cat(1,data_struct.shuff_label);
feature_cell = cat(2,num2cell(classpcolor),file_names,num2cell(shuff_label));

[sorted_cell,sorted_idx] = sortrows(feature_cell);
% allocate memory for a sorting vector
sorting_vector = zeros(length(data_struct),1);

% define the flag to prevent re-editting the structure
flag = 1;

% get the dataset name
stim_name = data_struct(1).name;
%scan for the p17b
if contains(stim_name, 'p17b')
    %if it's p17b
    stim_labels = {'R','G','B','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'R CK','UV CK','R GR','UV GR','R FL','UV FL'};
end
% define the color scheme depending on the stimulus type
if contains(data_struct(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
else
    color_scheme = distinguishable_colors(6);
end

% define the colors for the dataset
if data_struct(1).subsample == 1
    dataset_colors = [0 0 0;0 0 0];
    dataset_labels = {'RGCs','RAs'};
elseif data_struct(1).subsample == 4
    dataset_colors = paths.afOT_colors;
    dataset_labels = {'AF10','Tectum'};
else
    dataset_colors = paths.afOT_colors;
    dataset_labels = {'AF10','Tectum'};

end
% get the number of datasets
num_data = size(data_struct,2);
% get the number of stimuli
stim_num = size(data_struct(1).class{1}{1},1);
% Plot the classifier results
close all
% define the fontsize
fontsize = 18;
% allocate memory to store the matrices to plot later
conf_cell = cell(length(data_struct),1);
% for all the data sets
for datas = 1:num_data
    
    % if there are more datasets
    if num_data > 2
        counter = 2-mod(ceil(datas/2),2);
    else
        counter = datas;
    end
    
    % load the variables of interest
    class_cell = data_struct(datas).class;
    num_regions = data_struct(datas).num_regions;
    reg_label = {data_struct(datas).region{:,2}};
    empty_vec = data_struct(datas).empty_vec;
    % load the parameters
    name = data_struct(datas).name;
    classpcolor = data_struct(datas).classpcolor;
    subsample = data_struct(datas).subsample;
    loo = data_struct(datas).loo;
    shuff_label = data_struct(datas).shuff_label;
    repeat_number = data_struct(datas).repeat_number;
    bin_width = data_struct(datas).bin_width;
    redec_num = data_struct(datas).redec_num;
    portion = data_struct(datas).portion;
    cluster_flag = data_struct(datas).cluster_flag;
    if contains(data_struct(datas).name,{'Syn','syn'})
        data_struct(datas).figure_name = 'RGCs';
    else
        data_struct(datas).figure_name = 'Tectum';
    end
    % define the suffix to save images with
    suffix = strjoin({'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'portion',num2str(portion),...
        'cluster',num2str(cluster_flag),'.eps'},'_');
end
%% 1) Plot the performances side by side by region
% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 1 && num_data == 4
    close all
    fontsize = 12;
%     % set a counter for the x coordinate
%     x_counter = 1;
    % for all the data sets
    for datas = 1:2:num_data
        h = figure;
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        % get the number of regions
        reg_number = size(class_cell_shuff,1);
        % for all the regions
        for regs = 1:reg_number
            plot(regs,class_cell_real{regs}{2},'o','MarkerEdgeColor',[0 0 0],'MarkerSize',3);
            hold on
            % calculate the exact accuracy
            mean_acc = mean(class_cell_real{regs}{2});
            plot(regs,mean_acc,'o','MarkerFaceColor',[0 0 0],...
                'MarkerEdgeColor',[0 0 0],'MarkerSize',3)
        end

        % for all the regions
        for regs = 1:reg_number
            plot(regs,class_cell_shuff{regs}{2},'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',3);
            mean_shuf = mean(class_cell_shuff{regs}{2});
            plot(regs,mean_shuf,'o','MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',3)
        end
        % rescale of the max value is too low
        results_perregion = horzcat(class_cell_real{:});
        if max(vertcat(results_perregion{2,:})) < 0.5
            y_lim = [0 0.5];
        else
            y_lim = [0 1];
        end

        set(gca,'TickLength',[0 0],'LineWidth',2)
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8 8])
        set(gca,'XTick',1:regs,'XTickLabels',data_struct(datas).region(:,2),'FontSize',fontsize,...
            'XTickLabelRotation',45,'TickLabelInterpreter','none')
        ylabel('Accuracy','FontSize',fontsize)
        set(gca,'YLim',y_lim,'XLim',[0 reg_number + 1])
        
        % plot the shuffle line    
        % get the number of categories
        num_category = size(class_cell_real{1}{1},1);
        plot(get(gca,'XLim'),[1/num_category,1/num_category],'r--','LineWidth',1)
        box off

    end
end
%% 2) Plot the performances per stimulus as a matrix across data sets

if subsample == 2 && num_data == 4
     close all
     h = figure;
     % allocate memory for the performance matrix (stim by shuffle or not)
     performance_matrix = zeros(stim_num,num_data/2,2);
     % initialize data counter
     data_counter = 1;
    % for all the data sets
    for datas = 1:2:num_data
       
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;

        % calculate the mean confusion matrix
        mean_acc = mean(class_cell_real{1}{1},3);
        % get the performances for each stimulus
        stim_ave = diag(mean_acc)./sum(mean_acc,2);
        
        mean_shuf = mean(class_cell_shuff{1}{1},3);
        stim_shuffle = diag(mean_shuf)./sum(mean_shuf,2);
        
        % load them into the matrix
        performance_matrix(:,data_counter,1) = stim_ave;
        performance_matrix(:,data_counter,2) = stim_shuffle;
        % update the data counter
        data_counter = data_counter + 1;
        
    end
    
    % convert the rows to the stimulus color
    % allocate memory for the output matrix
    performance_colored = zeros(stim_num,num_data/2,4,2);
    %         performance_matrix = normr_1(performance_matrix,1);
    % for all the stimuli
    for stim = 1:stim_num
        if contains(stim_name, 'p17b')
            switch stim
                case 1
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(stim,:,:);
                case 2
                    
                    performance_colored(stim,:,1,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,2,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(stim,:,:);
                case 3
                    
                    performance_colored(stim,:,1,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(num_data/2,1,2);
                case 4
                    
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(num_data/2,1,2);
            end
            performance_colored(stim,:,4,:) = performance_matrix(stim,:,:);
        else
            switch stim
                case 1
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(1,:,:);
                    performance_colored(stim,:,4,:) = ones(num_data/2,1,2);
                case 2
                    
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,4,:) = performance_matrix(stim,:,:);
            end
        end
    end
    
    subplot(2,1,1)

    imagesc(performance_matrix(:,:,1))
    set(gca,'TickLength',[0 0],'XTick',[],'FontSize',fontsize)
    set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels)
    axis square
    set(gca,'CLim',[0 1])  
    
    subplot(2,1,2)
    imagesc(performance_matrix(:,:,2))
    set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels)
    set(gca,'TickLength',[0 0])
    set(gca,'CLim',[0 1])
    set(gca,'XTick',1:2,'XTickLabels',{data_struct([1 3]).figure_name},'FontSize',fontsize,...
        'XTickLabelRotation',45,'TickLabelInterpreter','none')
    axis square
    
end
%% 3) Plot classification over time

if num_data >= 4

    close all
    fontsize = 15;
       
    % allocate memory to store the values for stats
    mean_cell = cell(2,4);
    % set a counter for the x coordinate
    x_counter = 1;
    
    % define the data vector order depending on the protocol
    if contains(data_struct(1).name,'p17b')
        datas_vector = [3 1];
        f1 = figure;
        f2 = figure;
        f3 = figure;
        f4 = figure;
        
        fig_vector = [f1,f2,f3,f4];
    elseif contains(data_struct(1).name,'p8')
        datas_vector = [1:2:12];
        f1 = figure;
        f2 = figure;
        f3 = figure;
        f4 = figure;
        f5 = figure;
        f6 = figure;
        
        fig_vector = [f1,f2,f3,f4,f5,f6];
    end
        
    % for all the data sets
    for datas = datas_vector
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        
        % get the number of regions
        region_number = size(class_cell_real,1);
        % for all the regions
        for region = 1:region_number
            % reshape to put exp reps and class reps together
            real_pred = squeeze(mean(reshape(class_cell_real{region}{4}==class_cell_real{region}{5},...
                40,[],3,data_struct(datas).repeat_number),3));
            real_mean = mean(real_pred,3);

            real_sem = std(real_pred,0,3)./sqrt(size(real_pred,3));

            shuff_pred = squeeze(mean(reshape(class_cell_shuff{region}{4}==class_cell_shuff{region}{5},...
                40,[],3,data_struct(datas).repeat_number),3));
            shuff_mean = mean(shuff_pred,3);
            shuff_sem = std(shuff_pred,0,3)./sqrt(size(shuff_pred,3));
            
            % store the matrices for later calculations

            mean_cell{x_counter,1} = real_mean;
            mean_cell{x_counter,2} = shuff_mean;
            
            % also store the "pref" indexes
            mean_cell{x_counter,3} = (squeeze(sum(real_pred==1,1)-sum(real_pred~=1,1)))./...
                (squeeze(sum(real_pred==1,1)+sum(real_pred~=1,1)));

            % get the number of stimuli
            stim_num = size(real_mean,2);
            % for all the stimuli
            for stim = 1:stim_num
                if sum(ismember([3,5,7,9.11],datas))==1
                    color = dataset_colors(2,:);
                    trace = '-';
                else
                    color = dataset_colors(1,:);
                    trace = '-';
                end
                figure(fig_vector(stim));
                shadedErrorBar(1:size(real_mean,1),real_mean(:,stim),real_sem(:,stim),{'color',color,'linestyle',trace})
                hold on
                shadedErrorBar(1:size(shuff_mean,1),shuff_mean(:,stim),shuff_sem(:,stim))
                hold on
                set(gca,'YLim',[0 1.1])

                box off
                set(gca,'TickLength',[0 0])

            end
            set(gcf,'Name',data_struct(datas).region{region,2})

            % update the counter
            x_counter = x_counter + 1;
            
        end

    end

end
%% 4) Plot the success ratio

close all
% define the stimulus order
stim_order = [1 5 2 6 3 7 4 8];
% define the colors
stim_colors = [repmat(dataset_colors(2,:),4,1);repmat(dataset_colors(1,:),4,1)];

% get the real ratio
real_ratio = cat(1,mean_cell{:,3});

% plot
h = plotSpread(real_ratio(stim_order,:)','distributionColors',stim_colors(stim_order,:),'showMM',4);

% for all the colors
for colors = 1:length(stim_order)
    if mod(colors,2) == 0
        marker = 'o';
    else
        marker = 'o';
    end
    set(h{1}(colors),'markersize',5,'marker',marker)
end

set(h{3},'XTick',[])
set(h{2}(1),'Color','k','LineWidth',1)
set(h{2}(2),'Color','k','LineWidth',1)


% for all the pairs
for pair = 1:2:8
    ranksum(real_ratio(stim_order(pair),:),real_ratio(stim_order(pair+1),:))
end
%% 5) Plot varying numbers of neurons

% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 4 && num_data == 4
    close all
    fontsize = 18;
    % define the marker size
    marker_size = 2;
    
            
    % get the cell vector
    cell_vector = [5 10 20 40 80 100 150];
    
    % get the number of cell groups used
    cell_groups = length(cell_vector);
    
    % allocate memory to store the values for stats
    mean_cell = cell(2,cell_groups);
    % set a counter for the x coordinate
    x_counter = 1;
    h = figure;
    % for all the data sets
    for datas = [3 1]
        
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;

        
        % for all the cell groups
        for cells = 1:cell_groups
        
            % calculate the exact accuracy
            mean_acc = mean(class_cell_real{1}{2,cells});
            % store the accuracies for stats
            mean_cell{(datas+1)/2,cells} = class_cell_real{1}{2,cells};

            std_acc = std(class_cell_real{1}{2})./sqrt(size(class_cell_real{1}{2},1));
            hold on
            errorbar(cell_vector(cells),mean_acc,std_acc,'o','MarkerFaceColor',dataset_colors((datas+1)/2,:),...
                'MarkerEdgeColor',dataset_colors((datas+1)/2,:),'Color',dataset_colors((datas+1)/2,:),'MarkerSize',marker_size)
            mean_shuf = mean(class_cell_shuff{1}{2,cells});

            std_shuf = std(class_cell_shuff{1}{2})./sqrt(size(class_cell_real{1}{2},1));
            errorbar(cell_vector(cells),mean_shuf,std_shuf,'o',...
                'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
                'Color',[0.5 0.5 0.5],'MarkerSize',marker_size)

            % update the counter
            x_counter = x_counter + 1;
        end

    end
    % build the label vector
    cell_labels = string(cell_vector);
    cell_labels(end) = 'All';
    set(gca,'YLim',[0 1],'XTick',cell_vector(1:2:end),'XTickLabels',cell_labels(1:2:end),'XTickLabelRotation',45)

end
% Test the differences statistically

% for all the cell groups
for cells = 1:cell_groups
    p = ranksum(mean_cell{1,cells},mean_cell{2,cells})
end