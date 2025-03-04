function [groupstatesflag,index_tran,all_index,all_p2p_period,Op,Or]=distinguish_state_onlyP_succession(groupstates,positions,orientations,groupcentroid)    
numsteps= size(groupstates,1);
processing_skipvalue = 1;
groupstatesflag=zeros(numsteps,1);
for step =1:processing_skipvalue:numsteps
   positions_step=squeeze(positions(step,:,:))';
   orientations_step=squeeze(orientations(step,:,:));
   groupcentroid_step=squeeze(groupcentroid(step,:))';

   Op(step)=calculate_Op(orientations_step);
   Or(step)=calculate_Or(positions_step,groupcentroid_step,orientations_step);
   
   if all(groupstates(step,:)=='Polarized')
       groupstatesflag(step,1)=1;
   elseif all(groupstates(step,:)=='Milling  ')
       groupstatesflag(step,1)=2;
   elseif all(groupstates(step,:)=='Swarm    ')
       groupstatesflag(step,1)=3;
   elseif all(groupstates(step,:)=='-other   ')
       groupstatesflag(step,1)=4;
   end
%     if Op(step)>0.65
%         groupstatesflag(step,1)=1;
%     else
%         groupstatesflag(step,1)=4;
%     end
end
flag=groupstatesflag;
flag(flag~=4)=0;
flag(flag==4)=1;
dflag=[flag(1,:);flag(1:end-1,:)]-flag;
index_tran_0=find(dflag==-1);
index_tran_1=find(dflag==1)-1;
if index_tran_0(1)>index_tran_1(1)
    index_tran_1(1,:)=[];
end
% if index_tran_0(end)>index_tran_1(end)
%     index_tran_0(end,:)=[];
% end
if index_tran_0(end)>index_tran_1(end) & length(index_tran_0) > length(index_tran_1)
    index_tran_1(end+1,:)=size(groupstatesflag,1);
end
index_tran=[index_tran_0,index_tran_1]';

all_index1 = index_tran(1,:) - 1;
all_index2 = index_tran(2,:) + 1;
if all_index2(1)~=1
    all_index2 = [1 all_index2];
end
if all_index1(end)~=length(groupstatesflag)
    all_index1 = [all_index1 length(groupstatesflag)];
end
all_index = [all_index2;all_index1];
all_index = [all_index index_tran];
[~,index] = sort(all_index(1,:),'ascend');
all_index = all_index(:,index);

all_flag = groupstatesflag(all_index(2,:));

period_length = all_index(2,:)-all_index(1,:)+1;
% if sum(period_length)~=length(groupstatesflag)
%     error('Error!');
% end
% if sum(sum(index_tran - all_index(:,all_flag==4)))~=0
%     error('Error!');
% end

all_p2p_period = all_index(:,all_flag==1);
all_m2m_period = all_index(:,all_flag==2);

end