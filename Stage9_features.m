% Calculate calcium response standard features per color

% load the data
clearvars
close all
Paths
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;

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
%% 1) Average Gain plots
close all
% define the fontsize
fontsize = 20;
if contains(data(1).name,'p17b')
    % for all of the datasets
    for datas = 1:num_data

        % get the corresponding calcium data
        plot_matrix = data(datas).delta_norm;
        
        figure
        
        violins = violinplot(plot_matrix,[],'Bandwidth',0.1,'ViolinAlpha',0.5,'ShowData',true);
        % for all the stimuli
        for stim = 1:data(datas).stim_num
            violins(stim).ViolinColor = cone_color_scheme(stim,:);
            violins(stim).EdgeColor = cone_color_scheme(stim,:);
            violins(stim).MedianPlot.SizeData = 10;
            violins(stim).MedianPlot.MarkerEdgeColor = cone_color_scheme(stim,:);
        end
        set(gca,'XTick',[],'TickLength',[0 0])
        axis tight
        ylabel('Gain (a.u.)','FontName','Arial')

        [~,tbl,stats] = anova2(squeeze(plot_matrix),1,'off');
        s = multcompare(stats,'CType','bonferroni','Display','off');
    end
end
%% 2) Find the gain patterns

if contains(data(1).name,'p17b')
    
    close all
    % allocate memory to store the matrices
    type_cell = cell(num_data+1,3);
    
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
        [pattern_counts,sort_idx] = sort(pattern_counts,'descend');
        
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
        
        % store the matrix
        type_cell{datas,1} = pattern;
        type_cell{datas,2} = pattern_counts./sum(pattern_counts);
        type_cell{datas,3} = pattern_full;
        
        % eliminate the patterns with only 1 instance
        elim_vector = pattern_counts<2;
        pattern_counts = pattern_counts(~elim_vector);
        pattern_full = pattern_full(~elim_vector,:,:);
        
        figure
        set(gcf,'Color','w')
        subplot(2,1,2)
        image(permute(pattern_full,[2 1 3]))

        set(gca,'XLim',[-0.2 size(pattern_full,1)])
        set(gca,'YScale','linear','XTick',[],'Visible','off')

        subplot(2,1,1)
        BarPlotBreak(pattern_counts,pattern_counts(2)*1.8,pattern_counts(1)*0.9,'Line',0.6,2)
        set(gca,'YScale','linear','XTick',[],'Visible','off')
        

        axis tight
    end
    
end