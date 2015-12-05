This software provides a package provides several modules for joint analysis of image and shape collections

% Load the parameters
load('data\parameters.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attention: The rendering set in Para (i.e., Para.rHeights, 
% Para.rWidhts) are tested on Windows machines. If you are using
% linux or Mac, please call the following function to update 
% these two vectors:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Para = update_para(Para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module I: Pose estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module II: Pair-wise Image Shape Alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Module II Perform image-shape alignment for one image. It     
  selects a cloest shape, deform it to fit the image object. Note      
  that the pose estimation module is packed in the following main 
  function 

  [Camera_opt, Shape_opt] = demo_i2s_align(...
    Image,...  % Input image
    Shapes,... % Input shapes
    Para);     % Parameters used in alignment and please refer to 
                the function body for details

 Output argments:
  'Shape_opt' : The optimized shape, which aligns with the input  
                image object
  'Camera_opt': The associated optimized camera configuration



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module III: Mix-and-Match Image-Based Modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Mix-and-match image-based modeling stands for the task of    
  reconstructing the underlying 3D model of an image object by 
  assembling parts from a small collection of relevant shapes.   
  The input consists of one image and a few shapes, each of which 
  is given by a few disconnected components (e.g., parts). We  
  assume the shapes are consistently aligned in a world 
  coordinate system. We assume the pose of the image object is   
  given. The main function is given by

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Module IV: Joint Image-Shape Correspondence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Joint image-shape correspondence module estimates dense pixel-  
  wise correspondences among a collection of images and a   
  collection of shapes. We assume (i)the camera poses of the   
  images are pre-computed (i.e, using Module I), and (ii) The   
  shapes are consistently aligned in a world coordinate system. 
  The main function is given by

  PairMatches = demo_joint_i2s_corres_main(...
    Images,...       % Input images
    ImageCameras,... % The camera configurations associated with  
                       each input image
    Shapes,...       % The input shapes (which are aligned in a 
                       world coordinate system)
    Top,...          % (#Images x #Shapes) is a sparse matrix 
                       which specifies the pairs of image-shape 
                       for matching.
    Para);           %

  Here Top (#Images x #Shapes) is a sparse matrix which specifies  
  the pairs of image-shape for matching. For remaining function 
  arguments, please refer to the function body for details.

  Warning: If you throw in (>100) images and/or (>100) shapes, 
           the computation may take very long. So apply this 
           module wisely.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% External libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  The software package uses the following external libraries:

  Piotr's Image & Video Matlab Toolbox
  http://vision.ucsd.edu/~pdollar/toolbox/doc/
  http://research.microsoft.com/en-us/downloads/389109f6-b4e8-404c-84bf-239f7cbf4e3d/

  as well as the SIFT flow package
  http://people.csail.mit.edu/celiu/SIFTflow/

  HOG descriptor: functions 'imResample', 'hog'
  Edge Map: function 'edgesDetect'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


If you use the code, please cite the following paper:

Single-View Reconstruction via Joint Analysis of Image and Shape Collections, Qixing Huang, Hai Wang, and Vladlen Koltun. ACM Transactions on Graphics 34(4), 2015


