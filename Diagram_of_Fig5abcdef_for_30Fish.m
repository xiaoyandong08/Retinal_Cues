function Diagram_of_Fig5abcdef_for_30Fish


trials={'30-fish/0084/'};
basedir = '';
record_length = 100;

succession_type = 'onlyP';

fixed_length = 500;
blind_angle =  0/360*pi*2;

binocularangle = blind_angle/2;
nPhi=512;
phi=linspace(-pi,pi,2*nPhi+1);
phi=phi(phi>-pi+binocularangle & phi<pi-binocularangle);
plot_visual = 0;
anis_factor = 0;

for i = 1 : length(trials)

    if strcmp(succession_type,'onlyP')
        [all_tracks_filt{i},all_Frame_matrix{i},all_Op{i},all_Or{i},all_positions{i},all_lefteye{i},all_righteye{i},all_tailcoords{i},all_orientations{i},all_Frame_index{i},all_Frame_start_end{i}] = ...
            Generate_all_Frame_matrix_for_Fish(trials{i},basedir,record_length,succession_type,fixed_length);
    end
    for l = 13
        tic

        [~,~,~,all_visual_network,~,~,~,~]=...
            get_visual_occupy_and_theta(blind_angle,all_positions{i}{l},all_lefteye{i}{l},all_righteye{i}{l},all_tailcoords{i}{l},all_orientations{i}{l},plot_visual);

        [group_size{i}(l),~,ave_order{i}(l),diff_sign_retina_consensus{i}{l},...
            corr_retina_consensus{i}{l},ave_spatial_value{i}{l},ave_distance_value{i}{l}] ...
            = Calculate_immediate_couple_of_a_track_with_Smooth(all_Frame_matrix{i}{l},all_tracks_filt{i}{l},anis_factor);

 
        time = toc;
        disp([trials{i} ', ' num2str(l) ' of total ' num2str(length(all_Frame_matrix{i})) ', t = ' num2str(time) 's'])
    end
end


trial_index = 1;
threshold_ICNcorr = -0.4;
interval_ICN = [-0.4 -0.5 -0.6];
anis_factor = 0;

ICN_colors = [hex2rgb('FF0000');hex2rgb('00FF00');hex2rgb('0000FF');hex2rgb('999999')];

for kk = 13 
    traj_index = kk;

    Frame_start_end = all_Frame_start_end{1}(:,traj_index);
    load Frame_Image_30-fish-0084-traj013.mat

    Frame_matrix = all_Frame_matrix{trial_index}{traj_index};
    tracks_filt  = all_tracks_filt{trial_index}{traj_index};
    positions = all_positions{trial_index}{traj_index};
    lefteye  = all_lefteye{trial_index}{traj_index};
    righteye = all_righteye{trial_index}{traj_index};
    tailcoords  = all_tailcoords{trial_index}{traj_index};
    orientations = all_orientations{trial_index}{traj_index};
    Op = all_Op{trial_index}{traj_index};
    Or = all_Or{trial_index}{traj_index};
    Frame_index = all_Frame_index{trial_index}(:,traj_index);
    
    numfish = size(positions,2);
    agent_num = numfish;

    plot_visual = 0;
    [~,visual_theta,visual_dist,visual_network,V,~,~,~]=...
            get_visual_occupy_and_theta(blind_angle,all_positions{trial_index}{traj_index},all_lefteye{trial_index}{traj_index},all_righteye{trial_index}{traj_index},all_tailcoords{trial_index}{traj_index},all_orientations{trial_index}{traj_index},plot_visual);
    

    [retina_dist_ij_eachStep,retina_angle_ij_eachStep] = Calculate_Retina_dist_of_2frame(anis_factor,all_Frame_matrix{trial_index}{traj_index},all_tracks_filt{trial_index}{traj_index});
    [all_retina_eu_dist,all_retina_O_dist,all_order,all_Vretina_order,V_retina] = Calculate_core_of_ICN(threshold_ICNcorr,all_Frame_matrix{trial_index}{traj_index},all_tracks_filt{trial_index}{traj_index},corr_retina_consensus{trial_index}{traj_index},retina_dist_ij_eachStep,interval_ICN);
    
    
    all_data = {all_order all_retina_O_dist all_Vretina_order all_retina_eu_dist};

    figure;
    set(gcf,'position',[103 99 850 450])

    y_labels = {'ave order','$\left < R_{ij}(\tau=1)\right >_t$','ave order of ${\hat{\bf v}}^{\textrm{retina}}_{ij}(t)$','$\left < D_i^{\textrm{retina}} \right >_t$'};
    subpos = [0.62 0.60 0.15 0.35
              0.84 0.60 0.15 0.35
              0.62 0.09 0.15 0.35
              0.84 0.09 0.15 0.35];
    for i = 1 : length(all_data)
        subplot('Position',subpos(i,:))
        boxplot(all_data{i}','symbol', '','Width',0.8,'Labels',{['\Omega_{ij}<' num2str(interval_ICN(1))],['\Omega_{ij}<' num2str(interval_ICN(2))],['\Omega_{ij}<' num2str(interval_ICN(3))],'others'})
        %boxplot(all_data{i}','symbol', '','Width',0.8)
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),ICN_colors(4-j+1,:),'FaceAlpha',.4);
        end
        if i==1 || i == 3
            ytickformat('%.2f')
        end
        if i==2
            ax = gca;
            ax.YAxis.Exponent = -3;
        end
        set(gca,'FontSize',14);
        set(gca,'TickLabelInterpreter','tex');
        xlim([0.4 4.6])
        if i==1
            ylabel(y_labels{i})
        else
            ylabel(y_labels{i},'Interpreter','latex')
        end
    end



    ICN_corr = corr_retina_consensus{trial_index}{traj_index};
    ICN_corr(ICN_corr>threshold_ICNcorr) = 0;

    for step = 53%1:size(all_Frame_matrix{trial_index}{traj_index},2)

        for focal_idx = 22%1 : 10%agent_num

            ICN_neighbor = find(ICN_corr(focal_idx,:)<0);
            ICN_index6 = find(ICN_corr(focal_idx,:)<=interval_ICN(1));
            ICN_index7 = find(ICN_corr(focal_idx,:)<=interval_ICN(2));
            ICN_index8 = find(ICN_corr(focal_idx,:)<=interval_ICN(3));
            rest_index = setdiff(1:agent_num,[ICN_index6 ICN_index7 ICN_index8 focal_idx]);
            
            colors = zeros(agent_num,3);
            if length(ICN_index6)>0; colors(ICN_index6,:) = repmat(ICN_colors(1,:),length(ICN_index6),1);end
            if length(ICN_index7)>0; colors(ICN_index7,:) = repmat(ICN_colors(2,:),length(ICN_index7),1);end
            if length(ICN_index8)>0; colors(ICN_index8,:) = repmat(ICN_colors(3,:),length(ICN_index8),1);end
            if length(rest_index)>0; colors(rest_index,:) = repmat(ICN_colors(4,:),length(rest_index),1);end
            colors(focal_idx,:) = hex2rgb('000000');
            colorsAlpha = 0.3*ones(1,size(colors,1));
            colorsAlpha(rest_index) = 0.5;            
            
            subplot('position',[0.04 0.35 0.5 0.60])

            imshow(squeeze(frame_Image(step,:,:,:)))
            hold on
            fsegs = cat(3,squeeze(positions(step,:,:)),squeeze(lefteye(step,:,:)),squeeze(tailcoords(step,:,:)),squeeze(righteye(step,:,:)));
            centre_fish=squeeze(mean(fsegs,3));

            plot_fish_body=cat(3,squeeze(fsegs),squeeze(fsegs(:,:,1)));
            for i=1:numfish
                if i~=focal_idx
                    plot(squeeze(plot_fish_body(i,1,:))',squeeze(plot_fish_body(i,2,:))','color',colors(i,:),'linewidth',1);
                else
                    plot(squeeze(plot_fish_body(i,1,:))',squeeze(plot_fish_body(i,2,:))','color',colors(i,:),'linewidth',2);
                end
                hold on
                %text(centre_fish(i,1),centre_fish(i,2),num2str(i),'color','k','HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
                plot(squeeze(positions(1:step,i,1)),squeeze(positions(1:step,i,2)),'-','color',[colors(i,:) 0.3],'linewidth',1);
            end

            for jj = 1 : numfish
                see_id = jj;
                if focal_idx ~= see_id
                    temp_theta=visual_theta{step,focal_idx}{see_id};
                    viewpoint= mean([squeeze(lefteye(step,focal_idx,:)),squeeze(righteye(step,focal_idx,:))],2);
                    orient = orientations(step,focal_idx);
                    %quiver(viewpoint(1),viewpoint(2),cos(orientations(step,jj)),sin(orientations(step,jj)),100,'color','r','LineWidth',1)
                    dist=sqrt(sum((centre_fish-viewpoint').^2,2));
                    if temp_theta(2)-temp_theta(1)~=0
                        for i = 1 : size(temp_theta,2)
                            f = [1 2 3];
                            v = [viewpoint';[viewpoint(1,:)+dist(see_id)*cos(temp_theta(1,i)+orient) viewpoint(2,:)+dist(see_id)*sin(temp_theta(1,i)+orient)];...
                                [viewpoint(1,:)+dist(see_id)*cos(temp_theta(2,i)+orient) viewpoint(2,:)+dist(see_id)*sin(temp_theta(2,i)+orient)]];
                            patch('Faces',f,'Vertices',v,'EdgeColor','n','FaceColor',colors(see_id,:),'FaceAlpha',colorsAlpha(see_id));
                        end
                    end
                end
            end
            hold off
            set(gca,'xtick',[],'ytick',[])
            %xlim(xlims);ylim(ylims)
            axis equal
            LA=axis;
            text((LA(2)-LA(1))*0.05+LA(1),(LA(4)-LA(3))*0.93+LA(3),['Traj=' num2str(traj_index) ', Focal=' num2str(focal_idx)],'FontSize',16)
            title(['Frame = ' num2str(step)],'fontsize',14)
            
            axes('Position',[0.4 0.82 0.13 0.12])
            box on
            h = histfit(corr_retina_consensus{trial_index}{traj_index}(:),50,'kernel');
            xbar = h(1).XData;
            ybar = h(1).YData;
            xlines = h(2).XData;
            ylines = h(2).YData;
            bar(xbar, ybar/sum(ybar), 'BarWidth',1,'EdgeColor','n','FaceColor',hex2rgb('F6921E'),'FaceAlpha',0.7)
            hold on;
            plot(xlines,ylines/sum(ybar), 'Color',hex2rgb('F6921E'), 'LineWidth',2)
            xlabel('\Omega_{ij}')
            ylabel('P')
            xline(-0.4)
            xlim([-1 1])
            set(gca,'fontsize',10)

            if step>1
                subplot('position',[0.04 0.09 0.5 0.24])
                box on

                visual_distance_ij = squeeze(visual_dist(step,focal_idx,:));
                visual_theta_ij = visual_theta{step,focal_idx};
                visual_reach_neig = squeeze(visual_network(step,focal_idx,:))';

                diff_aa = cell2mat(V_retina(focal_idx,step-1));
                diff_aa_unit1 = normalized_vector(diff_aa)*0.5;

                for jj = 1 : numfish
                    if focal_idx ~= jj & visual_reach_neig(jj)~=0
                        for jj1 = 1 : size(visual_theta_ij{jj},2)
                            one_theta = visual_theta_ij{jj}(:,jj1);
                            if one_theta(2)-one_theta(1)>0

                                quiver(mean(one_theta)-diff_aa_unit1(jj,2),visual_distance_ij(jj),diff_aa_unit1(jj,2),diff_aa_unit1(jj,1),'color','k','LineWidth',1,'AutoScale', 'off', 'AutoScaleFactor', 0.2,'ShowArrowHead','off');
                                hold on
                                if diff_aa_unit1(jj,2)>0
                                    scatter(mean(one_theta),visual_distance_ij(jj),50,'Marker','<','MarkerEdgeColor','k','MarkerFaceColor','k')
                                else
                                    scatter(mean(one_theta),visual_distance_ij(jj),50,'Marker','>','MarkerEdgeColor','k','MarkerFaceColor','k')
                                end
                            end
                            f = [1 2 3 4];
                            v = [[one_theta(1) 0];[one_theta(1) visual_distance_ij(jj)];
                                [one_theta(2) visual_distance_ij(jj)];[one_theta(2) 0]];
                            patch('Faces',f,'Vertices',v,'EdgeColor','n','FaceColor',colors(jj,:),'FaceAlpha',0.5);


                        end
                    end
                end
                hold off
                set(gca,'xtick',[-pi -pi/2 0 pi/2 pi],'xticklabels',{'-\pi', '-\pi/2','0', '\pi/2', '\pi'})
                xline(0)
                xlim([-pi,pi])
                set(gca,'XDir','reverse','YTickLabel','','fontsize',14)
                ylabel('Distance');
                xlabel('Ego-22s view \phi_{22} of a frame shown in panel a')
            end
        end
        im(step) = getframe(gcf);
    end
    %saveas(gcf,['30-fish-0084-traj' num2str(traj_index,'%03d') '.jpg'])

end

end


