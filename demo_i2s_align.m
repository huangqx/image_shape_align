function [Camera_opt, Shape_opt] = demo_i2s_align(Image, Shapes, Para, verbose)

% Pre-compute hog descriptors of rendered images
Cameras = cam_camera_sampling(Para); % camera simulation
hogRender = cam_shape_hog_dess(Shapes, Cameras, Para, verbose);


% Initialize the camera parameter
Camera_init = cam_pose_est_single(Image, Shapes, Cameras, hogRender, Para, verbose);

% Perform the image shape alignment
[Camera_opt, Shape_opt] = i2s_pose_shape_refine(Image,...
    Shapes{Camera_init.closetShapeId},...
    Camera_init, Para);

