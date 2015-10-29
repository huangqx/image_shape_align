This software provides a package for joint image-shape analyasis and reconstruction

It uses two external libraries from 
Piotr's Image & Video Matlab Toolbox
http://vision.ucsd.edu/~pdollar/toolbox/doc/
http://research.microsoft.com/en-us/downloads/389109f6-b4e8-404c-84bf-239f7cbf4e3d/

as well as the SIFT flow package
http://people.csail.mit.edu/celiu/SIFTflow/

HOG descriptor: functions 'imResample', 'hog'
Edge Map: function 'edgesDetect'


% Load a dataset
load('data\car.mat');
% Each dataset consists of a few test images and a few shapes
% The shapes are consisttly oriented: the y-axis is the up-right orientation, and the underlying reflectional plane should be the yz-plane, the centers of the shapes are [0,0,0]
% Each image is given an estimated bounding box (does not need to be very accurate)

% Load the parameters
load('data\parameters.mat');

% Perform image-shape alignment for one image. It selects a cloest shape, deform it to fit the image object
[Camera_opt, Shape_opt] = demo_i2s_align(Images{1}, Shapes, Para);

Tips:
1. The most time-consuming part is computing the hog descriptors of the rendered images. These descriptors maybe pre-computed once.

2. Using more shapes would boost of the performance.

If you use the code, please cite the following paper:

Single-View Reconstruction via Joint Analysis of Image and Shape Collections, Qixing Huang, Hai Wang, and Vladlen Koltun. ACM Transactions on Graphics 34(4), 2015


