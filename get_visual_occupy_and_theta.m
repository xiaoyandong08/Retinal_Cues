function [visual_angular,theta,visual_dist,visual_network,V,V_id,GCC_size,visual_reach_percent]=get_visual_occupy_and_theta(blind_angle,positions,lefteye,righteye,tailcoords,orientations,plot_visual)

binocularangle = blind_angle/2;
nPhi=512;
phi=linspace(-pi,pi,2*nPhi+1);
phi=phi(phi>-pi+binocularangle & phi<pi-binocularangle);

[numsteps,numfish,~] = size(positions);
visual_angular=zeros(numsteps,numfish,numfish);
visual_dist=zeros(numsteps,numfish,numfish);
visual_network=zeros(numsteps,numfish,numfish);
metric_dist=zeros(numsteps,numfish,numfish);

% numsteps = 500;
for step=1:numsteps
    fsegs = cat(3,squeeze(positions(step,:,:)),squeeze(lefteye(step,:,:)),squeeze(tailcoords(step,:,:)),squeeze(righteye(step,:,:)));
    centre_fish=squeeze(mean(fsegs,3))';
    for i=1:numfish
        orient = orientations(step,i);
        eye= mean([squeeze(lefteye(step,i,:)),squeeze(righteye(step,i,:))],2);
        %             eye= squeeze(positions(step,i,:));
        focal_id=i;
%         if step ==28 & i==4
%             plot_visual =0;
%         else
%             plot_visual =0;
%         end
        [visual_angular(step,i,:),theta{step,i},visual_dist(step,i,:),metric_dist(step,i,:),V(step,i,:),V_id(step,i,:)] = get_visual_v3(focal_id,eye,orient,fsegs,centre_fish,phi,plot_visual);
    end
    visual_network(step,:,:) = sign(squeeze(visual_dist(step,:,:)));
    
    network = squeeze(visual_network(step,:,:));
    G = digraph(network);
    [bins,binsizes] = conncomp(G);
    GCC_size(step) = max(binsizes)/numfish;
    
    visual_reach_percent(step) = mean(sum(network,2)/(numfish-1));
end

end

function [visual_angular,theta,visual_dist,metric_dist,V,V_id]=get_visual_v3(focal_id,viewpoint,orient,segments,centre_fish,phi,plot_visual)

colors = [    0.0159    0.5512    0.7937
         0         0    1.0000
         0         0    0.0476
    0.9841         0    1.0000
         0    1.0000         0
    0.4286         0    0.5714
    0.4444    0.0157         0
    1.0000    0.5827         0
    1.0000    0.0079    0.4286
    0.1746    0.5039         0];

V=zeros(size(phi));
V_id=zeros(size(phi));
V_dist=ones(size(phi))*inf;
dist=sqrt(sum((centre_fish-viewpoint).^2,1));
% rng(0);
% colors=rand(size(segments,1),3);
id_in_vision = [];
for i=1:size(segments,1)
    if i~=focal_id
        segments_i=squeeze(segments(i,:,:));
        dist_i=dist(1,i);
        diffs_i=segments_i-viewpoint;
        diffs_i(:,vecnorm(diffs_i)==0)=[];
        theta_i=atan2(diffs_i(2,:),diffs_i(1,:));
        temp=theta_i-theta_i(1);
        theta_i(abs(temp)>pi)=theta_i(abs(temp)>pi)+sign(theta_i(1))*2*pi;
        theta_i=[min(theta_i),max(theta_i)]-orient;
        theta_i=wrapToPi(theta_i);
        if theta_i(2)>theta_i(1)
            V(phi>theta_i(1) & phi<theta_i(2))=1;
            V_id(phi>theta_i(1) & phi<theta_i(2) & V_dist>dist_i)=i;
            V_dist(phi>theta_i(1) & phi<theta_i(2) & V_dist>dist_i)=dist_i;
            if sum(phi>theta_i(1) & phi<theta_i(2))>0
                id_in_vision = [id_in_vision i];
            end
        else
            V(phi>theta_i(1) | phi<theta_i(2))=1;
            V_id(phi>theta_i(1) | phi<theta_i(2) & V_dist>dist_i)=i;
            V_dist(phi>theta_i(1) | phi<theta_i(2) & V_dist>dist_i)=dist_i;
            if sum(phi>theta_i(1) | phi<theta_i(2))>0
                id_in_vision = [id_in_vision i];
            end
        end
    end
end

for i=1:size(segments,1)
    if i==focal_id
        theta{1,i}=zeros(2,1);
    else
        temp=V_id;
        temp(temp~=i)=0;
        temp(temp==i)=1;
        dtemp=temp(1,2:end)-temp(1,1:end-1);
        flag0=find(dtemp==1)+1;
        flag1=find(dtemp==-1);
        
        %theta{1,i}=[phi(flag0);phi(flag1)];
        index=find(dtemp~=0);
        if isempty(flag0) & isempty(flag1)
            if temp(1)==1
                theta{1,i}=[phi(1);phi(end)];
            else
                theta{1,i}=zeros(2,1);
            end
        else
            if temp(1)==1
                index_0=[phi(1),phi(flag0)];
            else
                index_0=phi(flag0);
            end
            if temp(end)==1
                index_1=[phi(flag1),phi(end)];
            else
                index_1=phi(flag1);
            end
            theta{1,i}=[index_0;index_1];
            
            if dtemp(index(1))==-1
                if phi(1)==-pi
                    theta{1,i}=[theta{1,i}(1,2:end);[theta{1,i}(2,2:end-1),theta{1,i}(2,1)]];
                end
            end
        end
    end
end
if plot_visual == 1
    figure
    set(gcf,'position',[15 582 536 775])
    subplot('position',[0.1 0.36 0.8 0.6])
    for i=1:size(segments,1)
        if i~=focal_id
            temp_theta=theta{1,i};
            temp_theta=temp_theta(:);
            for j=1:numel(temp_theta)
                hold on;
                plot([viewpoint(1,:),viewpoint(1,:)+unique(V_dist(V_id==i))*cos(temp_theta(j)+orient)],[viewpoint(2,:),viewpoint(2,:)+unique(V_dist(V_id==i))*sin(temp_theta(j)+orient)],'color',colors(i,:));

            end
        end
    end
    temp=cat(3,squeeze(segments),squeeze(segments(:,:,1)));
    for i=1:size(segments,1)
        hold on;
        plot(squeeze(temp(i,1,:))',squeeze(temp(i,2,:))','color',colors(i,:),'linewidth',2);
%         f = [1 2 3 4];
%         v = squeeze(segments(i,:,:))';
%         patch('Faces',f,'Vertices',v,'EdgeColor',colors(i,:),'FaceColor',colors(i,:),'LineWidth',2);
        text(segments(i,1,1)-5,segments(i,2,1),num2str(i),'color','k','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    end
    title(num2str(id_in_vision))
    box on;
    axis equal;
    
    subplot('position',[0.1 0.05 0.8 0.25])
    hold on
    plot(phi,V,'-k','linewidth',2)
    %plot(phi,V_id,'-r','linewidth',1)
    ylim([0 1.05])
end
visual_angular=zeros(1,size(segments,1));
for i=1:size(segments,1)
    visual_angular(1,i)=numel(find(V_id==i))/size(phi,2);
end
temp = visual_angular*size(phi,2);
temp(temp~=0) = temp(temp~=0) - 1;
visual_angular = temp*(phi(2)-phi(1));
% aa = cell2mat(theta);
% aa(2,:) - aa(1,:) - visual_angular1;
%     [~,index]=sort(visual_angular,'descend');
%     visual_angular(1,index)=[1:size(segments,1)];

visual_dist = dist;
visual_dist(visual_angular==0) = 0;

metric_dist = zeros(1,length(dist));
metric_dist(id_in_vision) = dist(id_in_vision);
end


function [visual_angular,theta,visual_dist]=get_visual_v2(focal_id,viewpoint,orient,segments,centre_fish,phi,plot_visual)

colors = [    0.0159    0.5512    0.7937
         0         0    1.0000
         0         0    0.0476
    0.9841         0    1.0000
         0    1.0000         0
    0.4286         0    0.5714
    0.4444    0.0157         0
    1.0000    0.5827         0
    1.0000    0.0079    0.4286
    0.1746    0.5039         0];

V=zeros(size(phi));
V_id=zeros(size(phi));
V_dist=ones(size(phi))*inf;
dist=sqrt(sum((centre_fish-viewpoint).^2,1));
% rng(0);
% colors=rand(size(segments,1),3);
for i=1:size(segments,1)
    if i~=focal_id
        segments_i=squeeze(segments(i,:,:));
        dist_i=dist(1,i);
        diffs_i=segments_i-viewpoint;
        diffs_i(:,vecnorm(diffs_i)==0)=[];
        theta_i=atan2(diffs_i(2,:),diffs_i(1,:));
        temp=theta_i-theta_i(1);
        theta_i(abs(temp)>pi)=theta_i(abs(temp)>pi)+sign(theta_i(1))*2*pi;
        theta_i=[min(theta_i),max(theta_i)]-orient;
        theta_i=wrapToPi(theta_i);
        if theta_i(2)>theta_i(1)
            V(phi>theta_i(1) & phi<theta_i(2))=1;
            V_id(phi>theta_i(1) & phi<theta_i(2) & V_dist>dist_i)=i;
            V_dist(phi>theta_i(1) & phi<theta_i(2) & V_dist>dist_i)=dist_i;
        else
            V(phi>theta_i(1) | phi<theta_i(2))=1;
            V_id(phi>theta_i(1) | phi<theta_i(2) & V_dist>dist_i)=i;
            V_dist(phi>theta_i(1) | phi<theta_i(2) & V_dist>dist_i)=dist_i;
        end
    end
end

for i=1:size(segments,1)
    if i==focal_id
        theta{1,i}=[];
    else
        temp=V_id;
        temp(temp~=i)=0;
        temp(temp==i)=1;
        dtemp=temp-circshift(temp,1);
        theta{1,i}=[phi(dtemp==1);phi(dtemp==-1)];
        index=find(dtemp~=0);
        if isempty(index)
            if temp(1)==1
                theta{1,i}=[phi(1);phi(end)];
            else
                theta{1,i}=[];
            end
        else
            if dtemp(index(1))==-1
                %theta{1,i}=[theta{1,i}(1,:);[theta{1,i}(2,2:end),theta{1,i}(2,1)]];
                theta{1,i}=[theta{1,i}(1,:);phi(end)]; % xiao revised
                % 这里有问题
            end
        end
    end
end
if plot_visual == 1
    figure
    for i=1:size(segments,1)
        if i~=focal_id
            temp_theta=theta{1,i};
            temp_theta=temp_theta(:);
            for j=1:numel(temp_theta)
                hold on;
                plot([viewpoint(1,:),viewpoint(1,:)+unique(V_dist(V_id==i))*cos(temp_theta(j)+orient)],[viewpoint(2,:),viewpoint(2,:)+unique(V_dist(V_id==i))*sin(temp_theta(j)+orient)],'color',colors(i,:));
            end
        end
    end
    temp=cat(3,squeeze(segments),squeeze(segments(:,:,1)));
    for i=1:size(segments,1)
        hold on;
        plot(squeeze(temp(i,1,:))',squeeze(temp(i,2,:))','color',colors(i,:),'linewidth',2);
        text(segments(i,1,1)-5,segments(i,2,1),num2str(i),'color','k','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    end
    title(['Focal id = ' num2str(focal_id)])
    grid on;
    box on;
    axis equal;
end
visual_angular=zeros(1,size(segments,1));
for i=1:size(segments,1)
    visual_angular(1,i)=numel(find(V_id==i))/size(phi,2);
end
%     [~,index]=sort(visual_angular,'descend');
%     visual_angular(1,index)=[1:size(segments,1)];
visual_dist = dist;
visual_dist(visual_angular==0) = 0;
end