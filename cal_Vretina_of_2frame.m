function [Vretina] = cal_Vretina_of_2frame(xyz_pre,v_pre,xyz_now,v_now,focal_id)

diff_pos_pre = xyz_pre - xyz_pre(focal_id,:);
diff_angle_pre = atan2(diff_pos_pre(:,2),diff_pos_pre(:,1)) - atan2(v_pre(focal_id,2),v_pre(focal_id,1));
diff_angle_pre(focal_id) = [];
diff_angle_pre(diff_angle_pre>pi) = diff_angle_pre(diff_angle_pre>pi)-2*pi;
diff_angle_pre(diff_angle_pre<-pi) = 2*pi + diff_angle_pre(diff_angle_pre<-pi);

diff_pos_now = xyz_now - xyz_now(focal_id,:);
diff_angle_now = atan2(diff_pos_now(:,2),diff_pos_now(:,1)) - atan2(v_now(focal_id,2),v_now(focal_id,1));
diff_angle_now(focal_id) = [];
diff_angle_now(diff_angle_now>pi) = diff_angle_now(diff_angle_now>pi)-2*pi;
diff_angle_now(diff_angle_now<-pi) = 2*pi + diff_angle_now(diff_angle_now<-pi);

diff_angle = diff_angle_now - diff_angle_pre;
index  = find(abs(diff_angle)>pi);
if length(index)>0
    diff_angle(index)=-1*sign(diff_angle(index)).*(2*pi-abs(diff_angle(index)));
end
Vretina = sign(diff_angle);
end
