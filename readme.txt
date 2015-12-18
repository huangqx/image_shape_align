%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% License and attribution:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This package provides software modules for image and shape analysis. The code is released under the MIT license and can be used for any purpose with proper attribution. The code accompanies the following paper, which should be cited in publications that use the provided modules:

Single-View Reconstruction via Joint Analysis of Image and Shape Collections
Qixing Huang, Hai Wang, and Vladlen Koltun
ACM Transactions on Graphics 34(4), 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  load('data\parameters.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attention: The rendering setting in Para (i.e., Para.rHeights,
% Para.rWidhts) is tested on Windows machines. If you are using
% Linux or Mac, please call the following function to update
% these two vectors:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Para = update_para(Para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the images and shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  load('data\car.mat');
  Note that the shapes are assumed to be consistently oriented in
  a world coordinate system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module I: Pose estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  The pose estimation module predicts the camera pose for the underlying object with respect to the world coordinate system of the shapes. This is done by sampling the camera poses, rendering the shapes from the camera poses, and comparing real images to rendered images. The most time-consuming part is image rendering, which is done by precomputation.

  Cameras = cam_camera_sampling(Para); % camera simulation
  hogRender = cam_shape_hog_dess(...   %
                Shapes,...  % The aligned input shapes
                Cameras,... % The sampled camera poses
                Para,...    % Please refer to the function body
                verbose);   % verbose = 1 if you want to print it
                              out


  Camera estimation can be done in two ways: for a single image
  in isolation or for a collection of images

  Camera_init = cam_pose_est_single(...
                  Image,...     % The input image object
                  Shapes,...    % The aligned input shapes
                  Cameras,...   % The sampled camera poses
                  hogRender,... % The rendered hog descriptors
                  Para,...  % Please refer to the function body
                  verbose); % verbose = 1 if you want to print it
                              out

  Cameras_init = cam_pose_est_joint(...
                  inputImages,...  % The input images
                  Shapes,...       % The aligned input shapes
                  cameraSamples,...% The sampled camera poses
                  hogRender,...    % The rendered hog descriptors
                  Para,...   % Please refer to the function body
                  verbose)   % verbose = 1 if you want to print
                               it out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module II: Pairwise image-shape alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  This module performs image-shape alignment for one image. It selects the closest shape and deforms it to fit the image object. Note that the pose estimation module is packed in the following main function:

  [Camera_opt, Shape_opt] = demo_i2s_align(...
    Image,...  % Input image
    Shapes,... % Input shapes
    Para);     % Parameters used in alignment. Please refer to the
                function body for details

 Output argments:
  'Shape_opt' : The optimized shape, which aligns with the input
                image object
  'Camera_opt': The associated optimized camera configuration

Demo:
   >load('data\chair.mat');
   >load('data\parameters.mat');
   >[Camera_opt, Shape_opt] = demo_i2s_align(Images{1}, Shapes, Para);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module III: Assembly-based reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Assembly-based reconstruction aims to recover the underlying 3D model of an image object by assembling parts from a small collection of relevant shapes. The input consists of one image and a few shapes, each of which comprises a number of components. We assume that the shapes are consistently aligned in a world coordinate system. We assume that the pose of the image object is given. The main function is

  [Shape_opt, Camera_opt] = demo_ibm_mix_and_match(...
    Image,...             % The input image
    Camera_init,...       % The initial camera configuration for
                            all input shapes
    Shapes,...            % The input shapes
    Para);                % Please refer to the functional body
                            for details

  Output argments:
    'Shape_opt' : The reconstructed shape in the world coordinate
                  system
    'Camera_opt': The associated optimized camera configuration

Demo:
   >load('data\mix_and_match.mat');
   >load('data\parameters.mat');
   >[Shape_opt, Camera_opt] = demo_ibm_mix_and_match(Image, Camera_init, Shapes, Para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module IV: Joint image-shape correspondence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Joint image-shape correspondence module estimates dense pixelwise correspondences among a collection of images and a collection of shapes. We assume that (i) the camera poses of the images are pre-computed (e.g., using Module I), and (ii) the shapes are consistently aligned in a world coordinate system. The main function is

  PairMatches = demo_joint_i2s_corres_main(...
    Images,...       % Input images
    ImageCameras,... % The camera configurations associated with
                       each input image
    Shapes,...       % The input shapes (aligned in a
                       world coordinate system)
    Top,...          % (#Images x #Shapes) is a sparse matrix
                       that specifies image-shape pairs
                       for matching.
    Para);           %

  Warning: If you use >100 images and/or >100 shapes,
           the computation may take a very long time.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  The software package uses the following external libraries:

  Piotr's Image & Video Matlab Toolbox:
  http://vision.ucsd.edu/~pdollar/toolbox/doc/
  http://research.microsoft.com/en-us/downloads/389109f6-b4e8-404c-84bf-239f7cbf4e3d/

  The SIFT flow package:
  http://people.csail.mit.edu/celiu/SIFTflow/

  The OPENGM2 structure predication package:
  http://hci.iwr.uni-heidelberg.de/opengm2/?l0=library

  HOG descriptor: functions 'imResample', 'hog'

  Edge map: function 'edgesDetect'

  OPENGM2: the trws algorithm

References:

Navneet Dalal, Bill Triggs: Histograms of Oriented Gradients for Human Detection. CVPR 2005

Piotr Dollár, C. Lawrence Zitnick: Fast Edge Detection Using Structured Forests. IEEE Trans. Pattern Anal. Mach. Intell. 37(8): 1558-1570 (2015)

Ce Liu, Jenny Yuen, Antonio Torralba: SIFT Flow: Dense Correspondence across Scenes and Its Applications. IEEE Trans. Pattern Anal. Mach. Intell. 33(5): 978-994 (2011)

Andres, B. and Beier T. and Kappes, J.H. : OpenGM: A C++ Library for Discrete Graphical Models. http://arxiv.org/abs/1206.0111