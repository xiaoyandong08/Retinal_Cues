function Diagram_of_Fig4abc_for_group05


%%
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
% Plot_order_trajectory_for_One_Frame_onlyTraj(Frame_matrix,T,tracks_filt)

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

%% 
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
trace_length = 30;
if is_generate_video == 1
    figure('Color','w');
    set(gcf,'Position',[103 99 1000 426])
    %set(gcf,'Position',[103 99 1241 651])
    subplot('position',[0.16 0.03 0.8 0.1]);
    hold on
    text(0.3,1,'From ego-24,s view, neighbor-j with','FontSize',14)
    xlim([0 1])
    ylim([0 1])
    scatter(0.62,1,200,'ro','filled')
    text(0.64,0.91,'\Omega_{24,j}<-0.6','FontSize',14)
    scatter(0.75,1,200,'go','filled')
    text(0.77,0.91,'\Omega_{24,j}<-0.7','FontSize',14)
    scatter(0.88,1,200,'bo','filled')
    text(0.90,0.91,'\Omega_{24,j}<-0.8','FontSize',14)
    scatter(0.62,0.3,200,'o','MarkerEdgeColor',hex2rgb('A6A8AB'),'MarkerFaceColor',hex2rgb('A6A8AB'))
    text(0.64,0.21,'others (\Omega_{24,j}>=-0.6)','FontSize',14)
    
    annotation('textarrow',[0.3,0.35],[0.135,0.135],'String','{\bf V}_{24,j}^{retina}(t)   ','FontSize',14,'LineWidth',2);
    axis off
end

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

    %%%%%%%%%% 生成video
    if is_generate_video == 1

        subplot('position',[0.05 0.2 0.6 0.8]);
        imagesc(cube_img_1center);
        axis equal
        xticks = [0, size(cube_img_1center,2)]; 
        yticks = [0, size(cube_img_1center,1)]; 
        xlim([0,size(cube_img_1center,2)]);
        ylim([0,size(cube_img_1center,1)]);

        [~,x_index] = arrayfun(@(x) min(abs(x-x_angle_range)),xticks_angle,'UniformOutput',false);
        x_index = cell2mat(x_index);
        [~,y_index] = arrayfun(@(x) min(abs(x-y_angle_range)),yticks_angle,'UniformOutput',false);
        y_index = cell2mat(y_index);
        set(gca, 'XTick',x_pixel_range(x_index),'YTick',y_pixel_range(y_index)); 
        %set(gca, 'xticklabel', xticklabels_angle,'yticklabel', yticklabels_angle,'TickLabelInterpreter','latex'); 
        x_labelArray = [xticklabels_angle;compose('%d',x_pixel_range(x_index));];
        y_labelArray = [yticklabels_angle;compose('%d',y_pixel_range(y_index));];
        %set(gca, 'xticklabel', strtrim(sprintf('%s\\newline%s\n', x_labelArray{:})),'yticklabel', strtrim(sprintf('%s\\newline%s\n', y_labelArray{:}))); 
        set(gca, 'xticklabel', x_labelArray(1,:),'yticklabel', y_labelArray(1,:));
        set(gca,'FontSize',14)
        title(sprintf('Ego = %d, Frame = %d',focal_idx,t));
        % 
        hold on;
        scatter(pixel_x(ICN_index6,end),pixel_y(ICN_index6,end),300,'r.');
        scatter(pixel_x(ICN_index7,end),pixel_y(ICN_index7,end),300,'g.');
        scatter(pixel_x(ICN_index8,end),pixel_y(ICN_index8,end),300,'b.');
        scatter(pixel_x(rest_index,end),pixel_y(rest_index,end),300,'.','MarkerEdgeColor',hex2rgb('999999'));

        quiver(pixel_x(ICN_index6,end),pixel_y(ICN_index6,end),-diff_aa_unit1(ICN_index6,2),-diff_aa_unit1(ICN_index6,1),0.5,'color','r','LineWidth',1);
        quiver(pixel_x(ICN_index7,end),pixel_y(ICN_index7,end),-diff_aa_unit1(ICN_index7,2),-diff_aa_unit1(ICN_index7,1),0.3,'color','g','LineWidth',2);
        quiver(pixel_x(ICN_index8,end),pixel_y(ICN_index8,end),-diff_aa_unit1(ICN_index8,2),-diff_aa_unit1(ICN_index8,1),0.2,'color','b','LineWidth',3);
        quiver(pixel_x(rest_index,end),pixel_y(rest_index,end),-diff_aa_unit1(rest_index,2),-diff_aa_unit1(rest_index,1),0.2,'color',hex2rgb('999999'),'LineWidth',2);

        hold off

        if t<=size(Frame_matrix,2)-1
            subplot('position',[0.735 0.63 0.25 0.3]);
            plot([1:t],order_Vretina_ICN_6,'r-','LineWidth',2);hold on
            plot([1:t],order_Vretina_ICN_7,'g-','LineWidth',2)
            plot([1:t],order_Vretina_ICN_8,'b-','LineWidth',2)
            plot([1:t],order_Vretina_rest,'-','Color',hex2rgb('A6A8AB'),'LineWidth',2)
            ylabel('order of {\bf V}_{24,j}^{retina}(t)')
            set(gca,'Xticklabel','')
            set(gca,'FontSize',14)
            hold off
        end

        subplot('position',[0.735 0.27 0.25 0.3]);
        plot([1:t],retina_eu_dist_ICN_6,'r-','LineWidth',2);hold on
        plot([1:t],retina_eu_dist_ICN_7,'g-','LineWidth',2)
        plot([1:t],retina_eu_dist_ICN_8,'b-','LineWidth',2)
        plot([1:t],retina_eu_dist_rest,'-','Color',hex2rgb('A6A8AB'),'LineWidth',2)
        ylabel('D_{24}^{retina}(t)')
        xlabel('t (frame)')
        set(gca,'FontSize',14)
        hold off
        

        im(t) = getframe(gcf);
    end
    %%%%%%%%%% Fig.3a
    if is_generate_snapshot == 1
        if t == 50
            figure;
            set(gcf,'Position',[103 144 949 566])
            subplot('position',[0.02 0.37 0.6 0.6]);
            imagesc(cube_img_1center);
            axis equal
            xticks = [0, size(cube_img_1center,2)]; 
            yticks = [0, size(cube_img_1center,1)]; 
            xlim([0,size(cube_img_1center,2)]);
            ylim([0,size(cube_img_1center,1)]);

            [~,x_index] = arrayfun(@(x) min(abs(x-x_angle_range)),xticks_angle,'UniformOutput',false);
            x_index = cell2mat(x_index);
            [~,y_index] = arrayfun(@(x) min(abs(x-y_angle_range)),yticks_angle,'UniformOutput',false);
            y_index = cell2mat(y_index);
            set(gca, 'XTick',x_pixel_range(x_index),'YTick',y_pixel_range(y_index)); 
            %set(gca, 'xticklabel', xticklabels_angle,'yticklabel', yticklabels_angle,'TickLabelInterpreter','latex'); 
            x_labelArray = [xticklabels_angle;compose('%d',x_pixel_range(x_index));];
            y_labelArray = [yticklabels_angle;compose('%d',y_pixel_range(y_index));];
            set(gca, 'xticklabel', strtrim(sprintf('%s\\newline%s\n', x_labelArray{:})),'yticklabel', strtrim(sprintf('%s\\newline%s\n', y_labelArray{:}))); 
            set(gca,'FontSize',14)
            title(sprintf('focal=%d, step=%d, spatial=%.2f, corrThresh=%.2f',focal_idx,t,ind_ave_spatial_value(focal_idx),threshold_ICNcorr));
            % 
            hold on;
            scatter(pixel_x(ICN_index6,end),pixel_y(ICN_index6,end),300,'r.');
            scatter(pixel_x(ICN_index7,end),pixel_y(ICN_index7,end),300,'g.');
            scatter(pixel_x(ICN_index8,end),pixel_y(ICN_index8,end),300,'b.');
            scatter(pixel_x(rest_index,end),pixel_y(rest_index,end),300,'.','MarkerEdgeColor',hex2rgb('999999'));

            quiver(pixel_x(ICN_index6,end),pixel_y(ICN_index6,end),-diff_aa_unit1(ICN_index6,2),-diff_aa_unit1(ICN_index6,1),0.5,'color','r','LineWidth',1);
            quiver(pixel_x(ICN_index7,end),pixel_y(ICN_index7,end),-diff_aa_unit1(ICN_index7,2),-diff_aa_unit1(ICN_index7,1),0.3,'color','g','LineWidth',2);
            quiver(pixel_x(ICN_index8,end),pixel_y(ICN_index8,end),-diff_aa_unit1(ICN_index8,2),-diff_aa_unit1(ICN_index8,1),0.2,'color','b','LineWidth',3);
            quiver(pixel_x(rest_index,end),pixel_y(rest_index,end),-diff_aa_unit1(rest_index,2),-diff_aa_unit1(rest_index,1),0.2,'color',hex2rgb('999999'),'LineWidth',2);

            hold off

        end

        %%%%%%%%%% Fig.3bc
        if t == maxstep
            
            figure;
            set(gcf,'Position',[911   375   306   379])
            %subaxis(3,1,2,'SpacingVertical',0.03,'MarginLeft',.15,'MarginRight',.05);%subplot('position',[0.35 0.04 0.25 0.26]);
            subplot(3,1,2)
            plot([2:t],order_Vretina_ICN_6,'r-','LineWidth',2);hold on
            plot([2:t],order_Vretina_ICN_7,'g-','LineWidth',2)
            plot([2:t],order_Vretina_ICN_8,'b-','LineWidth',2)
            plot([2:t],order_Vretina_rest,'-','Color',hex2rgb('999999'),'LineWidth',2)
            ylabel('order of {\bf v}_{24,j}^{retina}(t)')
            set(gca,'Xticklabel','')
            set(gca,'FontSize',14)

            %subaxis(3,1,3,'SpacingVertical',0.03,'MarginLeft',.15,'MarginRight',.05);%subplot('position',[0.35 0.04 0.25 0.26]);
            subplot(3,1,3)
            plot([1:t],retina_eu_dist_ICN_6,'r-','LineWidth',2);hold on
            plot([1:t],retina_eu_dist_ICN_7,'g-','LineWidth',2)
            plot([1:t],retina_eu_dist_ICN_8,'b-','LineWidth',2)
            plot([1:t],retina_eu_dist_rest,'-','Color',hex2rgb('999999'),'LineWidth',2)
            ylabel('D_{24}^{retina}(t)')
            xlabel('t (frame)')
            set(gca,'FontSize',14)
        end
    end

end

end


function Plot_order_trajectory_for_One_Frame_onlyTraj(Frame_matrix,all_time,tracks_filt)

Color = jet(size(Frame_matrix,2));
figure('units','inches','position',[5 5.4306 5.0556 4.5694]);
box on

set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
view([-162 78])
% plot3(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2),tracks_filt(Frame_matrix(Frame_matrix(:)>0),3),tracks_filt(Frame_matrix(Frame_matrix(:)>0),4),'.')
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
hold on;box on
for i = 1 : size(Frame_matrix,2)
    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
 
    % if i==1 & length(unique(sum(logical(Frame_matrix),1)))==1
    %     text(xyz(:,1),xyz(:,2),xyz(:,3),num2str([1:size(Frame_matrix,1)]'))
    % end
    % if i == size(Frame_matrix,2)
    %     scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    % end
%     im(i) = getframe;
end
for i = 1 : size(Frame_matrix,2)
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    plot3(xyz(24,1),xyz(24,2),xyz(24,3),'.','MarkerSize',10,'color',hex2rgb('F74461'));
end
Id = Frame_matrix(:,end);
xyz=tracks_filt(Id,2:4);
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
title(['Ego (red) = ' num2str(24)])

% for i = 1 : size(Frame_matrix,1)
%     Id = Frame_matrix(i,find(Frame_matrix(i,:)>0));
%     if length(Id)>0
%         xyz=tracks_filt(Id(end),2:4);
%         scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
%     end
% end
set(gca,'fontsize',14)
view(-67,77)
if exist('im')
    a=VideoWriter('0.5_0.45','MPEG-4');
    a.FrameRate = 10;
    a.Quality = 100;
    open(a);
    writeVideo(a,im);
    close(a)
end

end
