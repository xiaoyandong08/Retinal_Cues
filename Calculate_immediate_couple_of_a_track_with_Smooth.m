function [group_size,ave_curvature,ave_order,diff_sign_retina_consensus,corr_retina_consensus,...
    ave_spatial_value,ave_distance_value,retina_dist_ij,retina_angle_ij] = Calculate_immediate_couple_of_a_track_with_Smooth(Frame_matrix,tracks_filt,anis_factor)

group_size = size(Frame_matrix,1);
[~,ave_curvature,ave_order] = cal_ave_order(Frame_matrix,tracks_filt);

[all_cos_value,all_std_dist,cos_value_time] = Calculate_cosij('slope',Frame_matrix,tracks_filt);
[retina_dist_ij,retina_angle_ij,spatial_value_ij,distance_ij] = Calculate_Retina_dist_of_2frame(anis_factor,Frame_matrix,tracks_filt);

for i = 1 : size(Frame_matrix,1)
    for j = 1 : size(Frame_matrix,1)
        if i~=j
            all_retina_dist{i,j} = squeeze(retina_dist_ij(1:end,i,j))';
            all_retina_angle{i,j} = squeeze(retina_angle_ij(1:end,i,j))';
            all_spatial_value{i,j} = squeeze(spatial_value_ij(1:end,i,j))';
            all_distance{i,j} = squeeze(distance_ij(1:end,i,j))';

            aa = all_retina_dist{i,j};
            bb = cos_value_time{i,j}(2:end);
            cc = all_distance{i,j};

            aa = aa(2:end) - aa(1:end-1);
            bb = bb(2:end) - bb(1:end-1);
            cc = cc(2:end) - cc(1:end-1);

            aa = smooth(aa,0.1,'loess')';
            bb = smooth(bb,0.1,'loess')';
            cc = smooth(cc,0.1,'loess')';
    
            corr_retina_consensus(i,j) = corr(aa',bb','type','Spearman');
            diff_sign_retina_consensus(i,j) = sum(sign(aa).*sign(bb)==-1)/length(aa);
            diff_sign_dist_consensus(i,j) = sum(sign(cc).*sign(bb)==-1)/length(cc);

            % if diff_sign_retina_consensus(i,j)<0.6
            %     figure;plot(sign(aa).*sign(bb),'.')
            % end
            % diff_sign_retina_consensus(i,j) = dtw(aa,bb);
            % diff_sign_dist_consensus(i,j) = dtw(cc,bb);
        end
    end
end
ave_spatial_value = cell2mat(cellfun(@mean,all_spatial_value,'UniformOutput',false));
ave_distance_value = cell2mat(cellfun(@mean,all_distance,'UniformOutput',false));
diff_sign_retina_consensus(1:size(diff_sign_retina_consensus,1)+1:size(diff_sign_retina_consensus,1)^2)=nan;
diff_sign_dist_consensus(1:size(diff_sign_dist_consensus,1)+1:size(diff_sign_dist_consensus,1)^2)=nan;
corr_retina_consensus(1:size(corr_retina_consensus,1)+1:size(corr_retina_consensus,1)^2)=nan;

end