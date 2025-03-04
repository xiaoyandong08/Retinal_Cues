function Or=calculate_Or(positions_step,groupcentroid_step,orientations_step)
    numfish=size(positions_step,2);
    u_step=[cos(orientations_step);sin(orientations_step);zeros(1,numfish)];
    r_step=[(positions_step-groupcentroid_step)./vecnorm(positions_step-groupcentroid_step);zeros(1,numfish)];
    Or=vecnorm(mean(cross(u_step,r_step),2));
end