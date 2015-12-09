function [Shape_opt, Camera_opt] = demo_ibm_mix_and_match(...
    Image, Camera_init, Shapes, Para)
% Assembly-based reconstruction aims to recover the underlying 3D model of 
% an image object by assembling parts from a small collection of relevant 
% shapes. The input consists of one image and a few shapes, each of which 
% comprises a number of components. We assume that the shapes are consistently
% aligned in a world coordinate system. We assume that the pose of the 
% image object is given.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
%   Image: the   Input image
%   Camera_init: The initial camera position computed using the pose
%                estimation module
%   Shapes:      The input shapes (which are consistently aligned in a
%                world coordinate system) 
%   Para:        The input parameters (see for details)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
%   Shape_opt: The reconstructed shape in the same world-coordinate system
%   Camera_opt: The optimized corresponding camera configuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize part-base correspondences between the input image and the input
% shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    % Opt to do image-shape alignment to bring the rendered image close 
    % to the image object (this is always recommended, albeit it takes
    % longer time.
    for i = 1:length(Shapes)
        [ShapeCameras{i}, Shapes{i}] = i2s_pose_shape_refine(...
            Image, Shapes{i}, Camera_init, Para);
    end
else
    % For efficiency, we can always use the one obtained from pose estimation
    for i = 1:length(Shapes)
        ShapeCameras{i} = Camera_init;
    end
end

% Compute part-wise dense correspondences
for i = 1:length(Shapes)
    partCorresStructs{i} = image2shape_dense_corres(Image,...
        Shapes{i}, ShapeCameras{i}, Para); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform part-based assembly to obtain an initial reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Shape_init, partCorres] = part_retrieval(Image,...
    Shapes, partCorresStructs, Para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine the part-based reconstruction, so that the refined reconstruction
% aligns with the input image object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Shape_opt, Camera_opt] = part_based_image2shape_align(Image, Camera_init,...
    Shape_init, partCorres, Para);
