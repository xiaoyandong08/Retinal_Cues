function [XvRotyz,rthetaphi,RT]=rotate_to_vhat_i_RT_stable(xyz,vhat_i,i)
% vhat_i is the row vector
xyz_i=xyz(i,:);  
% 
v_copy=zeros(numel(vhat_i(:,1)),2);
v_copy(:,1)=atan2(vhat_i(:, 3),sqrt(vhat_i(:,1).^2 +vhat_i(:,2).^2));  % elevation
v_copy(:,2)=atan2(vhat_i(:, 2), vhat_i(:, 1));  % azimuth
Xv_focal=xyz-xyz(i,:);

%XvRot:rotate and translate each position around the focal individual index

Rz=[cos(-v_copy(2)), -sin(-v_copy(2)),0;
    sin(-v_copy(2)), cos(-v_copy(2)),0;
    0,0,1];

Ry=[cos(-v_copy(1)), 0, -sin(-v_copy(1));
    0, 1, 0;
    sin(-v_copy(1)),0, cos(-v_copy(1))];

% 
g_hat = [0,0,-1]';  
frame_x_hat = vhat_i';  
frame_y_hat = -cross(g_hat, frame_x_hat);  
frame_z_hat = cross(frame_x_hat, frame_y_hat);   
local_axes = [frame_x_hat,frame_y_hat,frame_z_hat];
RT=Ry * Rz * local_axes;
 


XvRotyz = Ry * Rz * Xv_focal';  % OK
XvRotyz=XvRotyz';

[az,elev,r] = cart2sph(XvRotyz(:,1),XvRotyz(:,2),XvRotyz(:,3));
rthetaphi = [r elev az];




end
