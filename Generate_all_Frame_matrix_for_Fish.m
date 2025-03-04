function [tracks_filt,Frame_matrix,Op,Or,all_positions,all_lefteye,all_righteye,all_tailcoords,all_orientations,index_tran,Frame_start_end] = Generate_all_Frame_matrix_for_Fish(trials,basedir,record_length,succession_type,fixed_length)

trial=trials;
outname=[basedir,trial,'basicdata-v3.mat'];

load(outname);

if strcmp(succession_type,'onlyP') | strcmp(succession_type,'onlyT')
    [groupstatesflag,index_tran,all_index,all_p2p_period]=distinguish_state_onlyP_succession(groupstates,positions,orientations,groupcentroid);
elseif strcmp(succession_type,'Board')
    [groupstatesflag,index_tran]                         =distinguish_state_board(positions,groupstates,groupcentroid,record_length);
elseif strcmp(succession_type,'fixed')
    [index_tran]                         =distinguish_state_fixed_length(positions,groupstates,groupcentroid,fixed_length);
elseif strcmp(succession_type,'1state')
    [groupstatesflag,index_tran]=distinguish_state_1state(groupstates,positions,orientations,groupcentroid);

end

if strcmp(succession_type,'onlyP')
    colum_num=size(all_p2p_period,2);
elseif strcmp(succession_type,'onlyT') | strcmp(succession_type,'Board') | strcmp(succession_type,'fixed') | strcmp(succession_type,'1state')
    colum_num=size(index_tran,2);
end
count = 1;
for i=1:colum_num
    i;
    load(outname);
    
    if strcmp(succession_type,'onlyP')
        index = all_p2p_period(:,i);
    elseif strcmp(succession_type,'onlyT') | strcmp(succession_type,'Board')  | strcmp(succession_type,'fixed') | strcmp(succession_type,'1state')
        index = index_tran(:,i);
    end
    state_step = [index(1) index(end)];
    
    xlims = [20.8538620080479,1919.56118297756];
    ylims = [11.2523828031018,1026.41764230681];
    
    if state_step(2)-state_step(1)+1>record_length
        
        groupcentroid=groupcentroid(state_step(1):state_step(2),:);
        groupheading=groupheading(:,state_step(1):min(state_step(2),size(groupheading,2)));
        groupheadingxy=groupheadingxy(state_step(1):min(state_step(2),size(groupheadingxy,1)),:);
        groupstates=groupstates(state_step(1):state_step(2),:);
        lefteye=lefteye(state_step(1):state_step(2),:,:);
        orientations=orientations(state_step(1):state_step(2),:);
        positions=positions(state_step(1):state_step(2),:,:);
        righteye=righteye(state_step(1):state_step(2),:,:);
        rotcoords=rotcoords(state_step(1):state_step(2),:,:);
        tailcoords=tailcoords(state_step(1):state_step(2),:,:);
        
        all_positions{count} = positions;
        all_lefteye{count} = lefteye;
        all_righteye{count} = righteye;
        all_tailcoords{count} = tailcoords;
        all_orientations{count} = orientations;
        
%         xboard_dis = 100;
%         xboard = find(squeeze(positions(:,:,1)) - xlims(1) < xboard_dis | squeeze(xlims(2) -positions(:,:,1))  < xboard_dis);
%         yboard_dis = 50;
%         yboard = find(squeeze(positions(:,:,2)) - ylims(1) < yboard_dis | squeeze(ylims(2) -positions(:,:,2))  < yboard_dis);
%         
%         if (length(xboard)>1)|(length(yboard)>1)
%             reach_board{tnum}{index_i,index_j}(count_transition_state(index_i,index_j)) = 1;
%         else
%             reach_board{tnum}{index_i,index_j}(count_transition_state(index_i,index_j)) = 0;
%         end
        
        
        Frame_start_end(:,count) = [state_step(1);state_step(2)];
        [tracks_filt{count},Frame_matrix{count},Op{count},Or{count}] = for_one_frame(positions,orientations,groupcentroid);
        count = count + 1;
    end
end

end

function [tracks_filt,Frame_matrix,Op,Or] = for_one_frame(positions,orientations,groupcentroid)

all_vxyz = cat(3,cos(orientations),sin(orientations));
[numsteps,numfish,~] = size(positions);
Frame_matrix = reshape([1:1:numfish*numsteps],[numsteps,numfish])';
individual_id =  [1:1:numfish*numsteps]';
pxyz = reshape(positions,[numfish*numsteps,2]);
vxyz = reshape(all_vxyz,[numfish*numsteps,2]);
time = repmat([1:numsteps]',1,numfish);
time = time(:);
tracks_filt = [individual_id [pxyz zeros(numfish*numsteps,1)] time [vxyz zeros(numfish*numsteps,1)]];

for i = 1 : numsteps
    Op(i)=calculate_Op(orientations(i,:));
    Or(i)=calculate_Or(squeeze(positions(i,:,:))',groupcentroid(i,:)',orientations(i,:));
end



end