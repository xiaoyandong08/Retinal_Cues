function Op=calculate_Op(orientations_step)
    u_step=[cos(orientations_step);sin(orientations_step);];
    Op=vecnorm(mean(u_step,2));
end