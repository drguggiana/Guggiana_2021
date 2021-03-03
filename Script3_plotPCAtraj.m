%% 0) Clean up and load data

clearvars
close all
Paths
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).clusters_path;

data = load_clusters(cluster_path);
% Calculate PCA trajectories

close all
% define the stimulus set to use
stim_set = 1;

dataset_colors = paths.afOT_colors;

%define the stim labels based on the paradigm
%extract the actual file name
%scan for the p17b
if contains(data(1).name,'p17b')
    %if it's p17b
    stim_labels = {'R','G','B','UV'};
    %define the plot colors
    plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];
    plot_marker = {'o', 'o','o', 'o',};
    % define the CCA parameters
    var_threshold = 0.5;
    min_dim = 3;
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'R CK','UV CK','R GR','UV GR','R FL','UV FL'};
    %define the plot colors
    plot_col = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    plot_marker = {'o', 'o', 's', 's', 'd', 'd'};
    % define the CCA parameters
    var_threshold = 0.8;
    min_dim = 3;
end

% define the region set to use
region_set = 1;
% define which regions to include in the analysis
switch region_set
    case 1
        tectum_regions = {'R-TcN','R-TcP'};
        af_regions = {'AF10'};
        
    case 2
        tectum_regions = {'R-TcN','R-TcP','R-Cb','R-Hb','R-Pt','L-TcN','L-TcP','L-Cb','L-Hb','L-Pt'};
        af_regions = {'AF4','AF5','AF8','AF9','AF10'};
        
end
% get the number of datasets
num_data = size(data,2);
% get the time and stim nums
time_num = data(1).time_num;
stim_num = data(1).stim_num;
%% 1) Align pca trajectories from different fish with CCA and then plot
close all

% allocate memory to store the aligned PCs
cca_cell = cell(num_data,1);
% define the fontsize
fontsize = 15;
% define the time vector to take
time_vector = 21:60;
%for all the fish
for datas = 1:num_data
    
    %show the current fish
    fprintf(strcat('Current dataset:',num2str(datas),'\r\n'))
    
    % define the regions to be considered (depending on the stimulus protocol)
    if contains(data(datas).name,{'syn','Syn'})
        region_list = af_regions;
        if region_set == 1
            fig_name = 'AF10';
        else
            fig_name = data(datas).figure_name;
        end
    else
        region_list = tectum_regions;
        if region_set == 1
            fig_name = 'Tectum';
        else
            fig_name = data(datas).figure_name;
        end
    end
    
    % get the stimuli to use
    stim_vector = pick_stim_vector(stim_set,stim_num,data(datas).name);
    
    raw_trace = data(datas).conc_trace;
    
    resh_trace = reshape(raw_trace,size(raw_trace,1),time_num,stim_num);
    %clip the edges away
    resh_trace = resh_trace(:,time_vector,:);
    
    % get the number of animals and the animal info
    fish_ori_all = data(datas).fish_ori;
    num_animals = length(unique(fish_ori_all(:,1)));
    
    % exclude the anatomy if it's not present
    if isempty(data(datas).anatomy_info)
        anatomy_info = [];
    else
        anatomy_info = data(datas).anatomy_info(:,1);
    end
    [region_data,num_regions] = region_split(resh_trace,anatomy_info,data(datas).name,1,region_list);
    % for all the regions
    for region = 1:num_regions
        % get the region indexer
        region_idx = region_data{region,3};
        
        pca_mat = zeros(size(resh_trace,2),3,stim_num);
        
        % run the pcas per animal and store
        % initialize the skipped fish counter
        skip_count = 1;
        figure
        
        % allocate a structure to store the pca data
        pca_struct = struct([]);
        % split between animals
        % for all the animals
        for animals = 1:num_animals
            % get the traces for this animal
            temp_mat = reshape(resh_trace(fish_ori_all(:,1)==animals&region_idx==1,:,:),...
                [],stim_num*length(time_vector))';
            % run the pca
            [pca_struct(animals).coeff,pca_struct(animals).score,pca_struct(animals).latent] = ...
                pca(temp_mat);
            
        end
        
        % initialize the fish identifier
        target_fish = 0;
        % also a counter
        fish_counter = 1;
        % run the CCA and align the spaces
        while target_fish == 0
            % get the first animal's pca results
            var_1 = pca_struct(fish_counter).latent;
            traj_1  = pca_struct(fish_counter).score;
            
            %determine the 60% variance threshold
            dim_thres = find(cumsum(var_1./sum(var_1))>var_threshold,1,'first');
            % if it didn't pass the criteria
            if isempty(dim_thres) || dim_thres < min_dim
                % advance the counter
                fish_counter = fish_counter + 1;
                % if there are no fish that fulfill the criteria, error
                % out
                if fish_counter == num_animals
                    error('No or only one fish fulfill the criterion')
                else
                    fprintf(strcat('Skipped fish',num2str(skip_count),'\r\n'))
                    skip_count = skip_count + 1;
                    continue
                end
            else
                target_fish = fish_counter;
            end
        end
        
        %apply the threshold
        traj_1 = traj_1(:,1:dim_thres);
        
        % for all the animals minus the first one
        for animals = target_fish + 1:num_animals
            
            
            %get the corresponding pca results
            var_2 = pca_struct(animals).latent;
            traj_2  = pca_struct(animals).score;
            
            try
                traj_2 = traj_2(:,1:dim_thres);
            catch
                fprintf(strcat('Skipped fish',num2str(skip_count),'\r\n'))
                skip_count = skip_count + 1;
                continue
            end
            % align the current animal to the first one
            [pca_struct(target_fish).A,pca_struct(animals).B,...
                pca_struct(animals).r,~,~,stat] = canoncorr(traj_1,traj_2);
        end
        
        % plot the first animal
        pca_struct(target_fish).new_space = pca_struct(target_fish).score(:,1:dim_thres);
        % for all the animals
        for animals = target_fish+1:num_animals
            if isempty(pca_struct(animals).B) || (size(pca_struct(target_fish).A,1)>size(pca_struct(animals).B,1))
                continue
            end
            % transform the space
            pca_struct(animals).new_space = pca_struct(animals).score(:,1:dim_thres)*pca_struct(animals).B/...
                pca_struct(target_fish).A;
            
        end
        % eliminate the fish that don't fulfill the criterion

        % calculate the average trajectory in the common space
        all_animals = reshape(mean(cat(3,pca_struct.new_space),3),length(time_vector),stim_num,[]);
        % normalize the units for the space
        all_animals = normr_1(all_animals,1);
        
        % for all the stimuli
        for stim = stim_vector
            % plot the aligned trajectories
            options = struct([]);
            options(1).line=1;
            options(1).threeD=1;
            options(1).marker=plot_marker{stim};
            
            plot_trajectory(all_animals(:,stim,1:3),plot_col(stim,:),options)
        end
        
        % store the cca data (only for 1 region
        cca_cell{datas} = pca_struct;
        % get the handle of the current axis
        current_axis = gca;
        % set the line width of the plots
        view(2)
        
        set(gca,'YLim',[0 0.6],'XLim',[-0.1 1])
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        
    end
    
end
%% 2) Plot the single animal projections

close all
% set the colors and markers
cmap = [1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1];
% define the plotting options
options = struct([]);
options(1).line=1;
options(1).threeD=1;
options(1).marker=plot_marker{stim};
% define the labels
labels = {'Tectum','AF10'};
% for all the datasets
for datas = 1:num_data
    
    % get the data
    pca_struct = cca_cell{datas};
    current_stim = cat(3,pca_struct.score);
    current_stim = current_stim(:,1:3,:);

    % concatenate across stimuli
    all_trajectories = permute(reshape(current_stim,length(time_vector),stim_num,3,[]),[1 3 4 2]);
    
    % get the number of fish
    fish_num = size(all_trajectories,3);
    % for all the fish
    for fish = 1:fish_num
        figure
        % for all the stimuli
        for stim = 1:stim_num
            
            % plot it
            plot_trajectory(squeeze(all_trajectories(:,:,fish,stim)),cmap(stim,:),options)
        end
        view(2)
        set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',fontsize)
        set(gcf,'Color','w')
        axis equal
        title(strjoin({labels{datas},'Fish',num2str(fish)},' '),'FontSize',fontsize,'Interpreter','None')
        xlabel('PC 1','FontSize',fontsize)
        ylabel('PC 2','FontSize',fontsize)
        zlabel('PC 3','FontSize',fontsize)

    end
    
end
%% 3) Calculate the distance between components over time
close all

% allocate memory for the distances from the 2 datasets
distance_cell = cell(num_data,1);

% for all the datasets
for datas = 1:num_data
    
    % get the data
    pca_struct = cca_cell{datas};
    current_stim = cat(3,pca_struct.new_space);
    current_stim = current_stim(:,1:3,:);
    stim_matrix = permute(reshape(current_stim,length(time_vector),stim_num,3,[]),[1 3 4 2]);
    
    
    % get the number of fish
    fish_num = size(stim_matrix,3);
    % get the combination of stimuli
    stim_combo = nchoosek(1:stim_num,2);
    % get the number of combinations
    number_combos = length(stim_combo);

    % allocate memory for the resulting distances
    distances = zeros(number_combos,size(stim_matrix,1),fish_num);
    % for all the combos
    for combos = 1:number_combos
        % get the corresponding stimuli
        stim1 = stim_matrix(:,:,:,stim_combo(combos,1));
        stim2 = stim_matrix(:,:,:,stim_combo(combos,2));
        % allocate memory for the fish data
        fish_mat = zeros(fish_num,size(stim_matrix,1));
        % for all the fish
        for fish = 1:fish_num
            fish_mat(fish,:) = (vecnorm(stim1(:,:,fish)-stim2(:,:,fish),2,2));
        end
        % load the results matrix
        distances(combos,:,:) = fish_mat';

    end
    % save the distance matrix
    distance_cell{datas} = distances;
end
% allocate a matrix to plot
plot_matrix = zeros(stim_num);

% plot the results
for datas = 1:num_data
    % load and normalize the distance matrix
    distances = distance_cell{datas};

    % for all the combos
    for combos = 1:number_combos

        switch datas
            case 1
                x_coord = stim_combo(combos,1);
                y_coord = stim_combo(combos,2);
            case 2
                x_coord = stim_combo(combos,2);
                y_coord = stim_combo(combos,1);
        end
        
        % plot the distribution median in the matrix
        distance = squeeze(normr_1(distance_cell{datas}(combos,:,:),1));
        plot_matrix(x_coord,y_coord) = median(reshape(distance(:,:),[],1),1);
        
    end
end

% plot the matrix
imagesc(plot_matrix)

axis square
set(gca,'TickLength',[0 0],'XTick',1:stim_num,'YTick',1:stim_num)
set(gca,'XTickLabels',stim_labels,'YTickLabels',stim_labels,'XTickLabelRotation',45)

cba = colorbar;
ylabel(cba,'PC Normalized Distance (a.u.)')
cmap = magma;
cmap(1,:) = [1 1 1];
colormap(cmap)

% allocae memory for the test
test_medians = zeros(number_combos,1);
% for all combos
for combos = 1:number_combos
    
    distance1 = squeeze(normr_1(distance_cell{1}(combos,:,:),1));
    distance2 = squeeze(normr_1(distance_cell{2}(combos,:,:),1));
    % run the test, correcting for multiple comparisons
    test_medians(combos) = ranksum(distance1(:),distance2(:)).*number_combos;
    
end