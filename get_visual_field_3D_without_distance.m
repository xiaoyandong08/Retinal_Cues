function [r,XvRot] = get_visual_field_3D_without_distance(v,xyz,index)
%v_copy=[theta,phi]
v_copy=zeros(numel(v(:,1)),2);
v_copy(:,1)=atan2(v(:, 3),sqrt(v(:,1).^2 +v(:,2).^2));
v_copy(:,2)=atan2(v(:, 2), v(:, 1));
Xv=xyz-xyz(index,:);

%XvRot:rotate and translate each position around the focal individual index
XvRot=zeros(size(Xv));
XvRot(:,1)=Xv(:,1)*cos(-v_copy(index,2))-Xv(:,2)*sin(-v_copy(index,2));
XvRot(:,2)=Xv(:,1)*sin(-v_copy(index,2))+Xv(:,2)*cos(-v_copy(index,2));
XvRot(:,3)=Xv(:,3);

Xv=XvRot;
XvRot(:,1)=Xv(:,1)*cos(-v_copy(index,1))-Xv(:,3)*sin(-v_copy(index,1));
XvRot(:,2)=Xv(:,2);
XvRot(:,3)=Xv(:,1)*sin(-v_copy(index,1))+Xv(:,3)*cos(-v_copy(index,1));

% r=[r,theta,phi]
r=cart2sph(XvRot);


end

function r=cart2sph(x)
r=zeros(numel(x(:,1)),3);
r(:,1)=sqrt(x(:,1).^2 +x(:,2).^2+x(:,3).^2);
r(:,2)=atan2(x(:, 3),sqrt(x(:,1).^2 +x(:,2).^2));
r(:,3)=atan2(x(:, 2), x(:, 1));
end
