function [curvature,ave_curvature,ave_order] = cal_ave_order(Frame_matrix,tracks_filt)
order = nan(1,size(Frame_matrix,2));
moment = nan(1,size(Frame_matrix,2));
for i = 1 : size(Frame_matrix,2)
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    if length(Id)>0
        v_xyz = tracks_filt(Id,6:8);
        v_xyz = v_xyz./vecnorm(v_xyz')';
        order(i) = norm(sum(v_xyz,1))/size(v_xyz,1);
        
        p_xyz = tracks_filt(Id,2:4);
        mean_p_xyz = mean(p_xyz,1);
        p_xyz = (p_xyz - mean_p_xyz);
        p_xyz = p_xyz./vecnorm(p_xyz')';
        moment(i) = norm(sum(cross(p_xyz,v_xyz),1))/length(Id);

    end
end
ave_order = nanmean(order);

h = 1;
for id = 1 : size(Frame_matrix,1)
    id_xyz = tracks_filt(Frame_matrix(id,:),2:4);
    r = id_xyz';
    r1 = gradient(r)./h;
    r2 = gradient(r1)./h;
    r1 = r1';r2 = r2';
    v = cross(r1,r2,2);
    c = vecnorm(v');
    d = vecnorm(r1');
    k = (c./(d.^3));
    curvature(id,:) = k;
end
all_curvature = nanmean(curvature,1);
ave_curvature = nanmean(all_curvature);
end