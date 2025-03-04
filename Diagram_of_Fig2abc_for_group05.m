function Diagram_of_Fig2abc_for_group05


%% 取数据
load group_05.mat
% tracks = one_tracks_filt;
BirdIDs=unique(tracks(:,1));
T=unique(tracks(:,5));
tracks_filt = tracks;
Frame_matrix = zeros(length(BirdIDs),length(T));
for i = 1 : length(BirdIDs)
    Frame_matrix(i,:) = find(tracks(:,1)==BirdIDs(i));
    if sum(tracks(Frame_matrix(i,:),5)-T)~=0
        error('Error')
    end
end
Plot_order_trajectory_for_One_Frame_onlyTraj(Frame_matrix,T,tracks_filt)

for i = 1 : size(Frame_matrix,2)
    v = tracks(Frame_matrix(:,i),6:8);
    xyz = tracks(Frame_matrix(:,i),2:4);
    for j = 1 : size(Frame_matrix,1)
        r{j,i} = get_visual_field_3D_without_distance(v,xyz,j);
    end
end
for i = 1 : size(Frame_matrix,2)-1
    for j = 1 : size(Frame_matrix,1)
        aa_pre = cell2mat(r(j,i));
        aa_next = cell2mat(r(j,i+1));
        %diff_r{j,i} = aa_next-aa_pre;
        diff_r{j,i} = cell2mat(arrayfun(@(x) visual_difference(2*pi,aa_pre(x,[2 3]),aa_next(x,[2 3])),1:size(aa_pre,1),'UniformOutput',false)');
    end
end
diff_r(:,end+1) = diff_r(:,end);

%% 计算ICN及其他指标
threshold_ICNcorr = -0.6;
anis_factor = 0;



[group_size,ave_curvature,ave_order,diff_sign_retina_consensus,corr_retina_consensus,ave_spatial_value,ave_distance_value] = Calculate_immediate_couple_of_a_track(Frame_matrix,tracks_filt,anis_factor);
ind_ave_spatial_value = nanmean(ave_spatial_value,1);


[all_cos_value,all_std_dist,cos_value_time] = Calculate_cosij('slope',Frame_matrix,tracks_filt);
[retina_dist_ij,retina_angle_ij,spatial_value_ij,distance_ij] = Calculate_Retina_dist_of_2frame(anis_factor,Frame_matrix,tracks_filt);

ICN_corr = corr_retina_consensus;
ICN_corr(ICN_corr>threshold_ICNcorr) = 0;
percent_ICN = nansum(ICN_corr~=0,1)/size(ICN_corr,1);


image_dir=('UE5_for_Retina_reconstruction');
files = dir(image_dir);
maxstep=300;
minstep=1;

focal_idx=24;
neigh_idx = setdiff([1:70],focal_idx);
ICN_neighbor = find(ICN_corr(focal_idx,:)<0);
ICN_index6 = find(ICN_corr(focal_idx,:)<=-0.6);
ICN_index7 = find(ICN_corr(focal_idx,:)<=-0.7);
ICN_index8 = find(ICN_corr(focal_idx,:)<=-0.8);
rest_index = setdiff(1:70,[ICN_index6 ICN_index7 ICN_index8 focal_idx]);

pattern_id=sprintf('focal%d',BirdIDs(focal_idx));

xticks_angle = [pi pi/2 0 -pi/2 -pi];
yticks_angle = [pi/2 0 -pi/2];
xticklabels_angle = {'\pi','\pi/2','0','-\pi/2','-\pi'};
yticklabels_angle = {'\pi/2','0','-\pi/2'};
% xticklabels_angle = {'$\pi$','$\frac{\pi}{2}$','0','$-\frac{\pi}{2}$','$-\pi$'};
% yticklabels_angle = {'$\frac{\pi}{2}$','0','$-\frac{\pi}{2}$'};

BD=[min(tracks_filt(:,2))-1 max(tracks_filt(:,2))+1 ...
    min(tracks_filt(:,3))-1 max(tracks_filt(:,3))+1 ...
    min(tracks_filt(:,4))-1 max(tracks_filt(:,4))+1];

is_generate_snapshot = 1;
is_generate_video = 1 - is_generate_snapshot;


for t=1:maxstep
    vhat_i_now =  normalized_vector(tracks_filt(Frame_matrix(focal_idx,t),6:8));
    xyz_t_now = tracks_filt(Frame_matrix(:,t),2:4);
    u =  normalized_vector(tracks_filt(Frame_matrix(:,t),6:8));
    dist_to_focal = vecnorm((xyz_t_now-xyz_t_now(focal_idx,:))');

    [a,dist_sort_index] = sort(dist_to_focal,'ascend');

    aa = cell2mat(r(focal_idx,t));

    selected = aa([1:70],2:3);
    others = selected;
    for mm=1:size(selected,1)
        for nn=1:size(others,1)
            retina_eu_dist(mm,nn) = norm(visual_difference(2*pi,selected(mm,:),others(nn,:)));
        end
    end
    retina_eu_dist_ICN_6(t) = mean(mean(retina_eu_dist(ICN_index6,ICN_index6)));
    retina_eu_dist_ICN_7(t) = mean(mean(retina_eu_dist(ICN_index7,ICN_index7)));
    retina_eu_dist_ICN_8(t) = mean(mean(retina_eu_dist(ICN_index8,ICN_index8)));
    retina_eu_dist_rest(t) = mean(mean(retina_eu_dist(rest_index,rest_index)));

    if t<=size(Frame_matrix,2)-1
        retina_O_dist = squeeze(retina_dist_ij(t,focal_idx,:));
        retina_O_dist_ICN_6(t) = nanmean(retina_O_dist(ICN_index6));
        retina_O_dist_ICN_7(t) = nanmean(retina_O_dist(ICN_index7));
        retina_O_dist_ICN_8(t) = nanmean(retina_O_dist(ICN_index8));
        retina_O_dist_rest(t) = nanmean(retina_O_dist(rest_index));
    end

    order_ICN_6(t) = norm(sum(u([ICN_index6 focal_idx],:),1))/size(u([ICN_index6 focal_idx],:),1);
    order_ICN_7(t) = norm(sum(u([ICN_index7 focal_idx],:),1))/size(u([ICN_index7 focal_idx],:),1);
    order_ICN_8(t) = norm(sum(u([ICN_index8 focal_idx],:),1))/size(u([ICN_index8 focal_idx],:),1);
    order_rest(t) = norm(sum(u([rest_index focal_idx],:),1))/size(u([rest_index focal_idx],:),1);

    if t<=size(Frame_matrix,2)-1
        diff_aa = cell2mat(diff_r(focal_idx,t));
        diff_aa_unit1 = normalized_vector(diff_aa);
        order_Vretina_ICN_6(t) = norm(sum(diff_aa_unit1(ICN_index6,:),1))/size(diff_aa_unit1(ICN_index6,:),1);
        order_Vretina_ICN_7(t) = norm(sum(diff_aa_unit1(ICN_index7,:),1))/size(diff_aa_unit1(ICN_index7,:),1);
        order_Vretina_ICN_8(t) = norm(sum(diff_aa_unit1(ICN_index8,:),1))/size(diff_aa_unit1(ICN_index8,:),1);
        order_Vretina_rest(t) = norm(sum(diff_aa_unit1(rest_index,:),1))/size(diff_aa_unit1(rest_index,:),1);
    end

    %%%%%%%%%%%%%%%%%%%%
    pattern=sprintf('step%d.',t);
    filteredStruct = arrayfun(@(x) ((contains(x.name,pattern))&&(contains(x.name,pattern_id))), files, 'UniformOutput', false);
    filteredStruct_arr=cell2mat(filteredStruct);
    filteredStruct_val=files(logical(filteredStruct_arr));
    cube_img_path=sprintf('%s/%s',filteredStruct_val.folder,filteredStruct_val.name);
    cube_img=imread(cube_img_path);
    cube_img_1center=zeros(size(cube_img),'uint8');
    cube_img_1center(:,round(size(cube_img,2)/8*6+1):end,:)=cube_img(:,1:round(size(cube_img,2)/8*2),:);
    cube_img_1center(:,1:round(size(cube_img,2)/8*6),:)=cube_img(:,round(size(cube_img,2)/8*2+1):end,:);
    % cube_img_1center=cube_img_1center(:,end:-1:1,:);

    x_pixel_range = [1:1:size(cube_img_1center,2)];
    y_pixel_range = [1:1:size(cube_img_1center,1)];
    x_angle_range = linspace(pi,-pi,x_pixel_range(end));
    y_angle_range = linspace(pi/2,-pi/2,y_pixel_range(end));

    [~,x_index] = arrayfun(@(x) min(abs(x-x_angle_range)),aa(:,3),'UniformOutput',false);
    pixel_x(:,t) = cell2mat(x_index);
    [~,y_index] = arrayfun(@(x) min(abs(x-y_angle_range)),aa(:,2),'UniformOutput',false);
    pixel_y(:,t) = cell2mat(y_index);

    
    %%%%%%%%%% 生成Fig.2abc
    if is_generate_snapshot == 1    
        if t == maxstep
            figure;
            set(gcf,'Position',[103 144 920 566])
            subplot('position',[0.02 0.37 0.6 0.6]);
            imagesc(cube_img_1center);
            axis equal
            xticks = [0, size(cube_img_1center,2)]; % 定义刻度
            yticks = [0, size(cube_img_1center,1)]; % 定义刻度
            xlim([0,size(cube_img_1center,2)]);
            ylim([0,size(cube_img_1center,1)]);

            [~,x_index] = arrayfun(@(x) min(abs(x-x_angle_range)),xticks_angle,'UniformOutput',false);
            x_index = cell2mat(x_index);
            [~,y_index] = arrayfun(@(x) min(abs(x-y_angle_range)),yticks_angle,'UniformOutput',false);
            y_index = cell2mat(y_index);
            set(gca, 'XTick',x_pixel_range(x_index),'YTick',y_pixel_range(y_index)); % 设置刻度
            %set(gca, 'xticklabel', xticklabels_angle,'yticklabel', yticklabels_angle,'TickLabelInterpreter','latex'); % 设置刻度标签
            x_labelArray = [xticklabels_angle;compose('%d',x_pixel_range(x_index));];
            y_labelArray = [yticklabels_angle;compose('%d',y_pixel_range(y_index));];
            set(gca, 'xticklabel', strtrim(sprintf('%s\\newline%s\n', x_labelArray{:})),'yticklabel', strtrim(sprintf('%s\\newline%s\n', y_labelArray{:}))); % 设置刻度标签
            set(gca,'FontSize',14)
            title(sprintf('focal=%d, step=%d, spatial=%.2f, corrThresh=%.2f',focal_idx,t,ind_ave_spatial_value(focal_idx),threshold_ICNcorr));
            % 画ICN678的轨迹
            hold on;
            retina_focal = squeeze(retina_dist_ij(:,focal_idx,:));
            [~,max_index] = max(nanmean(retina_focal,1));
            max_index = 62;
            scatter(pixel_x(max_index,:),pixel_y(max_index,:),40,'.','MarkerEdgeColor',hex2rgb('F9EC31'));

            [~,min_index] = min(corr_retina_consensus(focal_idx,:));
            min_index = 18;
            scatter(pixel_x(min_index,:),pixel_y(min_index,:),40,'.','MarkerEdgeColor',hex2rgb('EC1C24'));

            hold off

            subplot('position',[0.70 0.735 0.25 0.18]);
            plot([1:t-1],retina_focal(:,max_index),'LineWidth',2,'Color',hex2rgb('262626'))
            ylabel(['R_{24,' num2str(max_index) '}(t,\tau=1)']);
            xlabel('t (frame)')
            set(gca,'FontSize',14,'YColor',hex2rgb('262626'))

            subplot('position',[0.70 0.43 0.25 0.18]);
            aa = squeeze(retina_dist_ij(1:end,focal_idx,min_index))';
            bb = cos_value_time{focal_idx,min_index}(2:end);

            aa = aa(2:end) - aa(1:end-1);
            bb = bb(2:end) - bb(1:end-1);
            corr(aa',bb','type','Spearman')

            yyaxis left
            plot(aa,'-','LineWidth',2)
            ylabel(sprintf('gradient of R_{%d,%d}',focal_idx,min_index))
            xlabel('t (frame)')
            yyaxis right
            plot(bb,'-','LineWidth',2)
            ylabel(sprintf('gradient of C_{%d,%d}',focal_idx,min_index))
            set(gca,'fontsize',14)
        end
    end

end

end


