function [all_retina_eu_dist,all_retina_O_dist,all_order,all_Vretina_order,diff_r] = Calculate_core_of_ICN(threshold_ICNcorr,Frame_matrix,tracks_filt,corr_retina_consensus,retina_dist_ij,interval_ICN)

agent_num = size(Frame_matrix,1);

for i = 1 : size(Frame_matrix,2)
    v = tracks_filt(Frame_matrix(:,i),6:8);
    xyz = tracks_filt(Frame_matrix(:,i),2:4);
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

for i = 1 : size(Frame_matrix,2)-1
    v_pre = tracks_filt(Frame_matrix(:,i),6:8);
    xyz_pre = tracks_filt(Frame_matrix(:,i),2:4);

    v_now = tracks_filt(Frame_matrix(:,i+1),6:8);
    xyz_now= tracks_filt(Frame_matrix(:,i+1),2:4);
    for j = 1: size(Frame_matrix,1)
        focal_id = j;
        [Vretina{j,i}] = cal_Vretina_of_2frame(xyz_pre,v_pre,xyz_now,v_now,focal_id);
        verify_diff(j,i) = sum(abs(Vretina{j,i} - sign(diff_r{j,i}(setdiff(1:agent_num,j),2))));
    end
end
if sum(abs(verify_diff(:)))>0
    error('Here is ERROR!!!!')
else
    disp('cal_Vretina_of_2frame function is Right!!!!')
end

ICN_corr = corr_retina_consensus;
ICN_corr(ICN_corr>threshold_ICNcorr) = 0;

for focal_idx = 1 : size(Frame_matrix,1)
    neigh_idx = setdiff([1:agent_num],focal_idx);

    ICN_neighbor = find(ICN_corr(focal_idx,:)<0);
    ICN_index6 = find(ICN_corr(focal_idx,:)<=interval_ICN(1));
    ICN_index7 = find(ICN_corr(focal_idx,:)<=interval_ICN(2));
    ICN_index8 = find(ICN_corr(focal_idx,:)<=interval_ICN(3));
    rest_index = setdiff(1:agent_num,[ICN_index6 ICN_index7 ICN_index8 focal_idx]);

    for t=1:size(Frame_matrix,2)
        vhat_i_now =  normalized_vector(tracks_filt(Frame_matrix(focal_idx,t),6:8));
        xyz_t_now = tracks_filt(Frame_matrix(:,t),2:4);
        u =  normalized_vector(tracks_filt(Frame_matrix(:,t),6:8));
        dist_to_focal = vecnorm((xyz_t_now-xyz_t_now(focal_idx,:))');

        aa = cell2mat(r(focal_idx,t));

        

        selected = aa([1:agent_num],2:3);
        others = selected;
        for mm=1:size(selected,1)
            for nn=1:size(others,1)
                retina_eu_dist(mm,nn) = norm(visual_difference(2*pi,selected(mm,:),others(nn,:)));
            end
        end
        retina_eu_dist(1:size(retina_eu_dist,1)+1:size(retina_eu_dist,1)^2) = nan;
        retina_eu_dist_ICN_6(t) = nanmean(nanmean(retina_eu_dist(ICN_index6,ICN_index6)));
        retina_eu_dist_ICN_7(t) = nanmean(nanmean(retina_eu_dist(ICN_index7,ICN_index7)));
        retina_eu_dist_ICN_8(t) = nanmean(nanmean(retina_eu_dist(ICN_index8,ICN_index8)));
        retina_eu_dist_rest(t) = nanmean(nanmean(retina_eu_dist(rest_index,rest_index)));
        clear retina_eu_dist
        
        if t<=size(Frame_matrix,2)-1
            retina_O_dist = squeeze(retina_dist_ij(t,focal_idx,:));
            retina_O_dist_ICN_6(t) = nanmean(retina_O_dist(ICN_index6));
            retina_O_dist_ICN_7(t) = nanmean(retina_O_dist(ICN_index7));
            retina_O_dist_ICN_8(t) = nanmean(retina_O_dist(ICN_index8));
            retina_O_dist_rest(t) = nanmean(retina_O_dist(rest_index));
        end
        
        % order_ICN_6(t) = norm(sum(u([ICN_index6 focal_idx],:),1))/size(u([ICN_index6 focal_idx],:),1);
        % order_ICN_7(t) = norm(sum(u([ICN_index7 focal_idx],:),1))/size(u([ICN_index7 focal_idx],:),1);
        % order_ICN_8(t) = norm(sum(u([ICN_index8 focal_idx],:),1))/size(u([ICN_index8 focal_idx],:),1);
        % order_rest(t) = norm(sum(u([rest_index focal_idx],:),1))/size(u([rest_index focal_idx],:),1);

        if t<=size(Frame_matrix,2)-1
            diff_aa = cell2mat(diff_r(focal_idx,t));
            diff_aa_unit1 = normalized_vector(diff_aa);
            
            % order_Vretina_ICN_6(t) = norm(sum(diff_aa_unit1(ICN_index6,:),1))/size(diff_aa_unit1(ICN_index6,:),1);
            % order_Vretina_ICN_7(t) = norm(sum(diff_aa_unit1(ICN_index7,:),1))/size(diff_aa_unit1(ICN_index7,:),1);
            % order_Vretina_ICN_8(t) = norm(sum(diff_aa_unit1(ICN_index8,:),1))/size(diff_aa_unit1(ICN_index8,:),1);
            % order_Vretina_rest(t) = norm(sum(diff_aa_unit1(rest_index,:),1))/size(diff_aa_unit1(rest_index,:),1);

            if length([ICN_index6])>=2
                order_Vretina_ICN_6(t) = norm(sum(diff_aa_unit1(ICN_index6,:),1))/size(diff_aa_unit1(ICN_index6,:),1);
            else
                order_Vretina_ICN_6(t) = nan;
            end
            if length([ICN_index7])>=2
                order_Vretina_ICN_7(t) = norm(sum(diff_aa_unit1(ICN_index7,:),1))/size(diff_aa_unit1(ICN_index7,:),1);
            else
                order_Vretina_ICN_7(t) = nan;
            end
            if length([ICN_index8])>=2
                order_Vretina_ICN_8(t) = norm(sum(diff_aa_unit1(ICN_index8,:),1))/size(diff_aa_unit1(ICN_index8,:),1);
            else
                order_Vretina_ICN_8(t) = nan;
            end
            if length([rest_index])>=2
                order_Vretina_rest(t) = norm(sum(diff_aa_unit1(rest_index,:),1))/size(diff_aa_unit1(rest_index,:),1);
            else
                order_Vretina_rest(t) = nan;
            end
        end

        if length([ICN_index6 focal_idx])>=2
            order_ICN_6(t) = norm(sum(u([ICN_index6 focal_idx],:),1))/size(u([ICN_index6 focal_idx],:),1);
        else
            order_ICN_6(t) = nan;
        end
        if length([ICN_index7 focal_idx])>=2
            order_ICN_7(t) = norm(sum(u([ICN_index7 focal_idx],:),1))/size(u([ICN_index7 focal_idx],:),1);
        else
            order_ICN_7(t) = nan;
        end
        if length([ICN_index8 focal_idx])>=2
            order_ICN_8(t) = norm(sum(u([ICN_index8 focal_idx],:),1))/size(u([ICN_index8 focal_idx],:),1);
        else
            order_ICN_8(t) = nan;
        end
        if length([rest_index focal_idx])>=2
            order_rest(t) = norm(sum(u([rest_index focal_idx],:),1))/size(u([rest_index focal_idx],:),1);
        else
            order_rest(t) = nan;
        end
    end

    mean_retina_eu_dist_ICN678{focal_idx}(1,1) = nanmean(retina_eu_dist_ICN_6);
    mean_retina_eu_dist_ICN678{focal_idx}(2,1) = nanmean(retina_eu_dist_ICN_7);
    mean_retina_eu_dist_ICN678{focal_idx}(3,1) = nanmean(retina_eu_dist_ICN_8);
    mean_retina_eu_dist_ICN678{focal_idx}(4,1) = nanmean(retina_eu_dist_rest);

    mean_retina_O_dist_ICN678{focal_idx}(1,1) = nanmean(retina_O_dist_ICN_6);
    mean_retina_O_dist_ICN678{focal_idx}(2,1) = nanmean(retina_O_dist_ICN_7);
    mean_retina_O_dist_ICN678{focal_idx}(3,1) = nanmean(retina_O_dist_ICN_8);
    mean_retina_O_dist_ICN678{focal_idx}(4,1) = nanmean(retina_O_dist_rest);

    mean_order_ICN678{focal_idx}(1,1) = nanmean(order_ICN_6);
    mean_order_ICN678{focal_idx}(2,1) = nanmean(order_ICN_7);
    mean_order_ICN678{focal_idx}(3,1) = nanmean(order_ICN_8);
    mean_order_ICN678{focal_idx}(4,1) = nanmean(order_rest);

    mean_order_Vretina_ICN678{focal_idx}(1,1) = mean(order_Vretina_ICN_6);
    mean_order_Vretina_ICN678{focal_idx}(2,1) = mean(order_Vretina_ICN_7);
    mean_order_Vretina_ICN678{focal_idx}(3,1) = mean(order_Vretina_ICN_8);
    mean_order_Vretina_ICN678{focal_idx}(4,1) = mean(order_Vretina_rest);

    clear retina_eu_dist_ICN_6 retina_eu_dist_ICN_7 retina_eu_dist_ICN_8 retina_eu_dist_rest
    clear retina_O_dist_ICN_6 retina_O_dist_ICN_7 retina_O_dist_ICN_8 retina_O_dist_rest
    clear order_ICN_6 order_ICN_7 order_ICN_8 order_rest 
    clear order_Vretina_ICN_6 order_Vretina_ICN_7 order_Vretina_ICN_8 order_Vretina_rest 

end
all_retina_eu_dist = cell2mat(mean_retina_eu_dist_ICN678);
all_retina_O_dist  = cell2mat(mean_retina_O_dist_ICN678);
all_order = cell2mat(mean_order_ICN678);
all_Vretina_order = cell2mat(mean_order_Vretina_ICN678);

end

