function Diagram_of_Fig2def_for_group05

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
% Plot_order_trajectory_for_One_Frame(Frame_matrix,T,tracks_filt,1)

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

%% 计算ICN及其他指标
threshold_ICNcorr = -0.6;
LF_correlation_threshold = 0.8;
anis_factor = 0;


[group_size,ave_curvature,ave_order,diff_sign_retina_consensus,corr_retina_consensus,...
    ave_spatial_value,ave_distance_value,retina_dist_ij,retina_angle_ij] = Calculate_immediate_couple_of_a_track(Frame_matrix,tracks_filt,anis_factor);
ind_ave_spatial_value = nanmean(ave_spatial_value,1);

%%%%%%%%%%%%%%%%%%%%%% Fig.2d %%%%%%%%%%%%%%%%%%%%%%
figure;
imagesc(corr_retina_consensus)
set(gca,'fontsize',12)
caxis([-1 0.6])
xlabel('Bird index')
ylabel('Bird index')
set(gca,'xtick',[10:10:70],'ytick',[10:10:70],'XTickLabelRotation',0)
set(gcf,'position',[835 455 253   197])
set(gca,'Position',[0.1449    0.1831    0.5913    0.7201])
title('CCD')
h = colorbar;
set(h,'Position',[0.76    0.1831    0.04    0.7201])

figure;
ind_corr_retina_consensus = nanmean(corr_retina_consensus,1);
[a,b] = sort(ind_corr_retina_consensus,'ascend');
stem(a,'MarkerSize',3,'Color',hex2rgb('0072BD'),'MarkerEdgeColor','n','MarkerFaceColor',hex2rgb('0072BD'))
xlabel('Bird index sorted by ascending \Omega_j')
ylabel('\Omega_j')
set(gca,'xticklabels','')
ylim([-0.8 0])
xlim([0 71])
set(gca,'fontsize',12)
set(gcf,'position',[835 455 253   87])

%%%%%%%%%%%%%%%%%%%%%% Fig.2e %%%%%%%%%%%%%%%%%%%%%%
aa = squeeze(nanmean(retina_dist_ij,1));
figure;
imagesc(aa)
set(gca,'fontsize',12)
xlabel('Bird index')
ylabel('Bird index')
set(gca,'xtick',[10:10:70],'XTickLabelRotation',0)
caxis([0 0.04])
set(gcf,'position',[835 455 253   197])
set(gca,'Position',[0.1449    0.1831    0.5913    0.7201])
title('$\left < R_{ij}(\tau=1)\right >$','Interpreter','latex')
h = colorbar;
set(h,'Position',[0.76    0.1831    0.04    0.7201])

figure;
ind_Rj = nanmean(squeeze(nanmean(retina_dist_ij,1)),1);
[a,b] = sort(ind_Rj,'ascend');
stem(a,'MarkerSize',3,'Color',hex2rgb('0072BD'),'MarkerEdgeColor','n','MarkerFaceColor',hex2rgb('0072BD'))
xlabel('Bird index sorted by ascending R_j')
ylabel('R_j')
set(gca,'xticklabels','')
% ylim([-0.8 0])
xlim([0 71])
set(gca,'fontsize',12)
set(gcf,'position',[835 455 265   87])

%%%%%%%%%%%%%%%%%%%%%% Fig.2f %%%%%%%%%%%%%%%%%%%%%%
aa = nanmean(retina_dist_ij,1);
aa = squeeze(aa);
aa = nanmean(aa,1);
bb = corr_retina_consensus;
bb = nanmean(bb,1);
figure;scatter(bb,aa,50,'filled');hold on;
x = bb;
y = aa;
[BS,xs]=func_cal_rlowess_bootstrap(x,y,10);
mean_BS=mean(BS,1);
Prctile_BS = prctile(BS,[3  97],1)';
h_area=area(xs',[Prctile_BS(:,1), (Prctile_BS(:,2)-Prctile_BS(:,1)) ],'FaceAlpha',0.5);
h_area(1).FaceColor = 'none';
h_area(2).FaceColor = hex2rgb('F7BECB');
h_area(1).EdgeColor = 'none';
h_area(2).EdgeColor = 'none';
h = plot(xs,mean_BS','-','linewidth',3,'color','r');
box on;
xlabel('\Omega_j');
ylabel('R_j')
set(gca,'fontsize',14,'xtick',[-0.7:0.1:-0.3])
xlim([-0.65 -0.35])
ylim([3*10^-3 10.5*10^-3])
set(gcf,'Position',[292   676   402   351])

end