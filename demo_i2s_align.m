function [Camera_opt, Shape_opt] = demo_i2s_align(Image, Shapes, Para)

% Initialize the camera parameter
Camera_init = cam_pose_est_single(Image, Shapes, Para, 1);

% Perform the image shape alignment
[Camera_opt, Shape_opt] = i2s_pose_shape_refine(Image,...
    Shapes{Camera_init.closetShapeId},...
    Camera_init, Para);

